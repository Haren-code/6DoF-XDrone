function [F_b, M_b] = dynamics(t, v_b, R_ib, omega_b, sim)
% DYNAMICS Compute total forces and moments for 6DoF rigid body
%
% Inputs:
%   v_b     - body-frame velocity [u; v; w]
%   R_ib    - rotation matrix inertial->body
%   params  - struct with .rho, .S, .m, .g, etc.
%   t       - current time (s)
%
% Outputs:
%   F_b     - total force in body frame
%   M_b     - total moment in body frame

    global SIM_DATA;

    % Compute velocity relative to wind (in body frame)
    if isfield(sim.aero,'V_wind_i')
        V_rel_b = v_b - R_ib * sim.aero.V_wind_i;
    else
        V_rel_b = v_b; % no wind
    end

    u_r = V_rel_b(1); 
    v_r = V_rel_b(2); 
    w_r = V_rel_b(3);
    V = norm(V_rel_b);

    if V < 1e-15
        alpha = 0; beta = 0;
    else
        alpha = atan2(w_r, u_r);   % AoA
        beta  = asin(v_r / V);     % sideslip
    end

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
    v_no_sdslip = V_no_sdslip_frame(2);
    w_no_sdslip = V_no_sdslip_frame(3);

    % Total angle of attack (in radians)
    if V < 1e-15
        alpha_equivalent = 0;
    else
        alpha_equivalent = atan2(w_no_sdslip, u_no_sdslip);
    end

    %%%%% Despining process
    % 1. Extract Roll (phi) robustly from R_ib
    % We use atan2(Ry, Rx) -> atan2(R(2,3), R(3,3))
    % This cancels out the cos(theta) term automatically.
    phi_despun = atan2(R_ib(2,3), R_ib(3,3));
    c_despun = cos(phi_despun);
    s_despun = sin(phi_despun);

    % 2. Create the "Un-roll" Matrix 
    % This is the Transpose of a standard X-axis rotation matrix.
    % It rotates Body Frame -> Non-Rotating Frame
    R_unroll = [1,   0,   0; 
                0,  c_despun,  s_despun; 
                0, -s_despun,  c_despun];

    % Dynamic pressure
    qbar = 0.5 * sim.aero.rho * V^2;

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

    % Gravity (in body frame)
    F_g = sim.prop.mass_total * [0; 0; sim.aero.g]; % inertial frame
    F_g_b = R_ib * F_g; % gravity in body frame

    if sim.options.propeller_on
        % --------- PROPELLER THRUST -------------
        % Prop RPM and relative rev/s
        n_rad_per_sec = abs(sim.prop.Prop.omega + omega_b(1));    % rad/s
        % Advance Ratio = axial velocity / (diameter*omega*2pi), cause it should be rev/s
        J = v_b(1) / (4*pi* n_rad_per_sec * sim.prop.Prop.radius);  
        % f_advance works like an efficiency - vector [1;0;0]
        % f_advance = max(0, 1 - params.Prop.k_V * J) * [1;0;0]; % reduces thrust at forward speed
        % f_advance = max(0, (1 - params.Prop.k_V  * J)^1.5)* [1;0;0];
        f_advance = max(0, 1 - sim.prop.Prop.k_V * J - sim.prop.Prop.k_V * J^2);

        % Thrust propeller
        T0_prop = sim.prop.Prop.k_T * n_rad_per_sec^2;
        T_prop = T0_prop * f_advance * [1;0;0]; 
    else
        T_prop = [0;0;0]; %xzylo
    end

    % Forces in body frame including gravity, without external forces
    F_b_no_external = F_aero_body + F_g_b + T_prop;

    % Add the F_external if it exist
    F_b = F_b_no_external + sim.aero.Fext(t);%[25*(1-u_r/sim.initial.V_mag);0;0];    

    % --------------------- Aerodynamic moments ---------------------------
    CoP_frac = sim.aero.Xzylo.CoP_frac(alpha_equivalent)/100;
    
    C_l = 0; % roll moment negligible for axisymmetric
    C_m = (sim.prop.percentage_CoG_total - CoP_frac) * ( C_L*c_alpha_eq + C_D*s_alpha_eq );
    C_n = 0; % 0 because we look in the no-sideslip frame

    M_no_sdslip = qbar * sim.prop.Area * [sim.prop.span * C_l
                                          sim.prop.chord * C_m
                                          sim.prop.span * C_n];

    % Rotate the moments around x axis to go back to body frame
    M_aero = R_no_sdslip' * M_no_sdslip;

    % --------- PROPELLER Reaction TORQUE -------------
    if sim.options.propeller_on
        Q_0 = sim.prop.Prop.k_Q * n_rad_per_sec^2;
        Q_mag = f_advance * Q_0;
    else
        Q_mag = [0; 0; 0];
    end

    T_friction = sim.prop.f_coeff*(pi*sim.aero.rho*omega_b(1)^2*sim.prop.span^4*sim.prop.chord)*[1;0;0];

    % Add to total forces
    M_b_no_external = M_aero + Q_mag - T_friction;

    % External moments addition
    M_b = M_b_no_external + sim.aero.Mext(t);

    % ------------------ LOGGING DATA ------------------
    SIM_DATA.t(end+1)  = t;
    SIM_DATA.V(end+1)  = V;
    SIM_DATA.V_no_sdslip(end+1) = norm(V_no_sdslip_frame); % should be equal to V from above
    SIM_DATA.u_no_sdslip(end+1) = u_no_sdslip;
    SIM_DATA.v_no_sdslip(end+1) = v_no_sdslip; % should be zero
    SIM_DATA.w_no_sdslip(end+1) = w_no_sdslip;

    SIM_DATA.Lift_sdslip(end+1) = -F_wind_frame(3); % I took lift and drag as positive values
    SIM_DATA.Drag_sdslip(end+1) = -F_wind_frame(1);
    SIM_DATA.Normal_force(end+1) = -F_aero_no_sdslip(3);
    SIM_DATA.Axial_force(end+1) = -F_aero_no_sdslip(1);

    SIM_DATA.Fx_aero_b(end+1) = F_aero_body(1); 
    SIM_DATA.Fy_aero_b(end+1) = F_aero_body(2);
    SIM_DATA.Fz_aero_b(end+1) = F_aero_body(3);

    SIM_DATA.Fx_b(end+1) = F_b(1); 
    SIM_DATA.Fy_b(end+1) = F_b(2);
    SIM_DATA.Fz_b(end+1) = F_b(3);

    F_nr = R_unroll' * F_b;
    SIM_DATA.Fx_nr(end+1) = F_nr(1);
    SIM_DATA.Fy_nr(end+1) = F_nr(2);
    SIM_DATA.Fz_nr(end+1) = F_nr(3);
    
    SIM_DATA.Mx_no_sdslip(end+1) = M_no_sdslip(1); % Should be zero 
    SIM_DATA.My_no_sdslip(end+1) = M_no_sdslip(2);
    SIM_DATA.Mz_no_sdslip(end+1) = M_no_sdslip(3); % Should be zero

    SIM_DATA.Mx_b(end+1) = M_b(1);
    SIM_DATA.My_b(end+1) = M_b(2);
    SIM_DATA.Mz_b(end+1) = M_b(3); 

    M_nr = R_unroll' * M_b;    % Moments
    SIM_DATA.Mx_nr(end+1) = M_nr(1);
    SIM_DATA.My_nr(end+1) = M_nr(2);
    SIM_DATA.Mz_nr(end+1) = M_nr(3);

    SIM_DATA.CL(end+1) = C_L; 
    SIM_DATA.CD(end+1) = C_D;

    SIM_DATA.Cm(end+1) = C_m;

    SIM_DATA.alpha(end+1) = alpha;
    SIM_DATA.beta(end+1)  = beta;
    SIM_DATA.no_sdslip_angle(end+1) = no_sdslip_angle;
    SIM_DATA.alpha_equivalent(end+1) = alpha_equivalent;

    SIM_DATA.CoP_fraction(end+1) = CoP_frac;

end