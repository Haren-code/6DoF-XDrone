function [F_b, M_b] = dynamics(t, v_b, omega_b, R_ib, sim)
% DYNAMICS Compute total forces and moments for 6DoF rigid body
%
% Inputs:
%   v_b      - body-frame velocity [u; v; w]
%   omega_b  - body-frame angular velocities [p; q; r]
%   R_ib     - rotation matrix inertial->body
%   sim      - struct with all info
%   t        - current time (s)
%
% Outputs:
%   F_b      - total force in body frame
%   M_b      - total moment in body frame


    %% ================= WIND AND VELOCITIES =================
    % Compute velocity relative to wind (in body frame)
    V_rel_b = v_b - R_ib * sim.prop.V_wind_i;

    % u_r = V_rel_b(1); %not used
    v_r = V_rel_b(2); 
    w_r = V_rel_b(3);
    V = norm(V_rel_b);

    %% ================= ANGLES & ROTATIONS =================
    %Alpha and beta angles are here, but actually never used in axisymetric
    %bodies, so it should be deleted, I don't think it helps with anything
    % if V < 1e-15
    %     alpha = 0; beta = 0;
    % else
    %     alpha = atan2(w_r, u_r);   % AoA
    %     beta  = asin(v_r / V);     % sideslip
    % end

    % Compute angle to rotate about body x so lateral component becomes zero
    if abs(w_r) < 1e-15 && abs(v_r) < 1e-15
        no_sdslip_angle = 0;
    else
        no_sdslip_angle = atan2(v_r, w_r); 
    end
    
    % Rotation matrix for sdslip about x-axis
    c_sdslip = cos(no_sdslip_angle);
    s_sdslip = sin(no_sdslip_angle);

    R_no_sdslip = [1 0 0
                   0 c_sdslip -s_sdslip
                   0 s_sdslip c_sdslip];

    % Express velocity in non-sideslip frame (should have v ≈ 0)
    V_no_sdslip_frame = R_no_sdslip * V_rel_b;
    u_no_sdslip = V_no_sdslip_frame(1);
    % v_no_sdslip = V_no_sdslip_frame(2); %not used
    w_no_sdslip = V_no_sdslip_frame(3);

    % Total angle of attack (in radians)
    if V < 1e-15
        alpha_equivalent = 0;
    else
        alpha_equivalent = atan2(w_no_sdslip, u_no_sdslip);
    end

    %% ================= AERODYNAMICS =================
    % Dynamic pressure
    qbar = 0.5 * sim.prop.rho * V^2;

    % Aerodynamic coefficients (based on total AoA only for axisymmetric)
    C_D = sim.aero.Xzylo.C_D(alpha_equivalent);
    C_L = sim.aero.Xzylo.C_L(alpha_equivalent);
    C_Y = 0; % MAYBE USE THIS FOR MAGNUS EFFECT LATER ON, NOT IMPLEMENTED YET

    % Forces in no-sideslip rotated frame (convention: drag opposes x_th, lift opposes z_th)
    F_wind_frame = qbar * sim.prop.Area * [-C_D; -C_Y; -C_L];

    % Forces must be rotated by alpha_equivalent to go in no-sideslip frame
    c_alpha_eq = cos(alpha_equivalent);
    s_alpha_eq = sin(alpha_equivalent);

    R_alpha_total = [c_alpha_eq 0 -s_alpha_eq
                     0 1 0
                     s_alpha_eq 0 c_alpha_eq];

    % Transformation from lift/drag to normal/axial
    F_aero_no_sdslip = R_alpha_total * F_wind_frame;
        
    % Rotate forces back to body frame
    F_aero_body = R_no_sdslip' * F_aero_no_sdslip;

    %% ================= PROPELLER =================
    if sim.options.propeller_on
        n_rad_per_sec = abs(sim.prop.Prop.omega + p);
        J = u ./ (4 * pi .* n_rad_per_sec .* sim.prop.Prop.radius);
        f_advance = max(0, 1 - sim.prop.Prop.k_V .* J - sim.prop.Prop.k_V .* J.^2);
        
        T_prop = f_advance * (sim.prop.Prop.k_T * n_rad_per_sec.^2) * [1;0;0];
        Q_prop = f_advance * (sim.prop.Prop.k_Q * n_rad_per_sec.^2) * [1;0;0];
    else
        T_prop = [0;0;0]; 
        Q_prop = [0;0;0];
    end

    %% ================= GRAVITY =================
    % Gravity (in body frame)
    F_g = sim.prop.mass_total * [0; 0; sim.prop.g]; % inertial frame
    F_g_b = R_ib * F_g; % gravity in body frame
    
    
    %% ================= TOTAL FORCES =================
    % Forces in body frame including gravity, external forces and propeller
    F_b = F_aero_body + F_g_b + T_prop + sim.aero.Fext(t);


    %% ================= TOTAL MOMENTS =================
    CoP_frac = sim.aero.Xzylo.CoP_frac(alpha_equivalent)/100;
    
    C_l = 0; % roll moment negligible for axisymmetric
    C_m = (sim.prop.percentage_CoG_total - CoP_frac) * ( C_L*c_alpha_eq + C_D*s_alpha_eq );
    C_n = 0; % 0 because we look in the no-sideslip frame

    M_no_sdslip = qbar * sim.prop.Area * [sim.prop.span  * C_l
                                          sim.prop.chord * C_m
                                          sim.prop.span  * C_n];

    % Rotate the moments around x axis to go back to body frame
    M_aero = R_no_sdslip' * M_no_sdslip;

    % Rotating drag / Friction torque
    T_friction = sim.prop.f_coeff*(pi*sim.prop.rho*omega_b(1)^2*sim.prop.span^4*sim.prop.chord)*[1;0;0];

    % Add all the moments
    M_b = M_aero + Q_prop - T_friction + sim.aero.Mext(t);

end