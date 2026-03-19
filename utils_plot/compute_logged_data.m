function SIM_DATA = compute_logged_data(t, x, sim)

    N = length(t);
    
    % Extract states into 1 x N arrays
    u = x(:, 1)'; v = x(:, 2)'; w = x(:, 3)'; 
    p = x(:, 4)'; % q_ang = x(:, 5); r = x(:, 6); 

    % Quaternions to Rotation Matrices (N x 4 input -> 3 x 3 x N output)
    quat = x(:, 7:10);
    q_norm = quatNorm(quat); 
    R_bi = quat2rotm(q_norm); 
    R_ib = pagetranspose(R_bi);
    
    %% ================= WIND AND VELOCITIES =================
    % V_rel_b = v_b - R_ib * sim.prop.V_wind_i;

    %First the R_ib * sim.prop.V_wind_i part
    V_wind_b_3D = pagemtimes(R_ib, sim.prop.V_wind_i); 
    V_wind_b = squeeze(V_wind_b_3D); % Flattens to 3 x N
    
    % Next part of the equation v_b - R_ib * sim.prop.V_wind_i;
    u_r = u - V_wind_b(1, :); % Needed for V norm
    v_r = v - V_wind_b(2, :);
    w_r = w - V_wind_b(3, :);

    V = sqrt(u_r.^2 + v_r.^2 + w_r.^2);
    
    %% ================= ANGLES & ROTATIONS =================
    no_sdslip_angle = zeros(1, N);
    valid_sdslip = abs(w_r) >= 1e-15 | abs(v_r) >= 1e-15;
    % no_sideslip_angle = arctan (v/w)
    no_sdslip_angle(valid_sdslip) = atan2(v_r(valid_sdslip), w_r(valid_sdslip));
    
    c_sd = cos(no_sdslip_angle); s_sd = sin(no_sdslip_angle);
    
    % V_no_sdslip_frame = R_no_sdslip * V_rel_b;
    u_no_sdslip = u_r;
    v_no_sdslip = v_r .* c_sd - w_r .* s_sd;
    w_no_sdslip = v_r .* s_sd + w_r .* c_sd;
    
    alpha_eq = zeros(1, N);
    valid_V = V >= 1e-15; 
    % alpha_equivalent = arctan(w_no_sdslip/u_no_sdslip);
    alpha_eq(valid_V) = atan2(w_no_sdslip(valid_V), u_no_sdslip(valid_V));
    
    %% ================= AERODYNAMICS =================
    qbar = 0.5 * sim.prop.rho * V.^2;
    
    C_D = sim.aero.Xzylo.C_D(alpha_eq); 
    C_L = sim.aero.Xzylo.C_L(alpha_eq);
    C_Y = zeros(1, N); 
    
    % F_wind_frame = qbar * sim.prop.Area * [-C_D; -C_Y; -C_L];
    F_wind_x = qbar .* sim.prop.Area .* (-C_D);
    F_wind_y = qbar .* sim.prop.Area .* (-C_Y);
    F_wind_z = qbar .* sim.prop.Area .* (-C_L);
    
    c_aeq = cos(alpha_eq); s_aeq = sin(alpha_eq);
    
    % F_aero_no_sdslip = R_alpha_total * F_wind_frame;
    F_ans_x = F_wind_x .* c_aeq - F_wind_z .* s_aeq;
    F_ans_y = F_wind_y;
    F_ans_z = F_wind_x .* s_aeq + F_wind_z .* c_aeq;
    
    % F_aero_body = R_no_sdslip' * F_aero_no_sdslip;
    Fx_aero_b =  F_ans_x;
    Fy_aero_b =  F_ans_y .* c_sd + F_ans_z .* s_sd;
    Fz_aero_b = -F_ans_y .* s_sd + F_ans_z .* c_sd;
    
    %% ================= GRAVITY =================
    % F_g_b = R_ib * F_g;
    F_g_z = sim.prop.mass_total * sim.prop.g;
    Fx_g_b = squeeze(R_ib(1,3,:))' .* F_g_z;
    Fy_g_b = squeeze(R_ib(2,3,:))' .* F_g_z;
    Fz_g_b = squeeze(R_ib(3,3,:))' .* F_g_z;
    
    %% ================= PROPELLER =================
    if sim.options.propeller_on
        R = sim.prop.Prop.radius;
    
        % Angular speed [rad/s]
        n_rad_per_sec = abs(sim.prop.Prop.omega + p);
    
        % Advance ratio (correct, vector-safe form)
        J = (pi * u_no_sdslip) ./ (n_rad_per_sec .* R);
    
        % Advance correction (your empirical model)
        f_advance = max(0, 1 ...
            - (J / 0.955).^1.5);   % thrust law
        
        f_advance_Q = max(0, 1 ...
            - 0.7 * (J / 0.955).^1.33);  % torque law
    
        % Precomputed dimensional constants
        kT_dim = sim.prop.Prop.k_T_0 * sim.prop.rho * (2*R)^4 / (2*pi)^2;
        kQ_dim = sim.prop.Prop.k_Q_0 * sim.prop.rho * (2*R)^5 / (2*pi)^2;
    
        % Final thrust & torque (vector form preserved)
        Tx_prop = kT_dim .* (n_rad_per_sec.^2) .* f_advance;
        Qx_prop = kQ_dim .* (n_rad_per_sec.^2) .* f_advance_Q;
    
    else
        Tx_prop = zeros(1, N);
        Qx_prop = zeros(1, N);
    end

    %% ================= EXTERNAL FORCES =================
    F_ext = zeros(3, N); M_ext = zeros(3, N);
    for i = 1:N
        F_ext(:, i) = sim.aero.Fext(t(i));
        M_ext(:, i) = sim.aero.Mext(t(i));
    end

    %% ================= TOTAL FORCES =================
    Fx_b = Fx_aero_b + Fx_g_b + Tx_prop + F_ext(1, :);
    Fy_b = Fy_aero_b + Fy_g_b + F_ext(2, :);
    Fz_b = Fz_aero_b + Fz_g_b + F_ext(3, :);
    
    %% ================= TOTAL MOMENTS =================
    CoP_frac = sim.aero.Xzylo.CoP_frac(alpha_eq) / 100;
    
    C_l = zeros(1, N);
    C_m = (sim.prop.percentage_CoG_total - CoP_frac) .* (C_L .* c_aeq + C_D .* s_aeq); 
    C_n = zeros(1, N);
    
    Mx_no_sd = qbar .* sim.prop.Area .* sim.prop.span  .* C_l;
    My_no_sd = qbar .* sim.prop.Area .* sim.prop.chord .* C_m;
    Mz_no_sd = qbar .* sim.prop.Area .* sim.prop.span  .* C_n;
    
    Mx_aero =  Mx_no_sd;
    My_aero =  My_no_sd .* c_sd + Mz_no_sd .* s_sd;
    Mz_aero = -My_no_sd .* s_sd + Mz_no_sd .* c_sd;
    
    T_friction = sim.prop.f_coeff .* (pi .* sim.prop.rho .* (p.^2) .* sim.prop.span^4 .* sim.prop.chord); 
    
    Mx_b = Mx_aero + Qx_prop - T_friction + M_ext(1, :);
    My_b = My_aero + M_ext(2, :);
    Mz_b = Mz_aero + M_ext(3, :);
    
    %% ================= PACK SIM_DATA =================
    % (Transposing everything back to Nx1 columns for plots)
    SIM_DATA.t = t; 
    SIM_DATA.V = V';
    SIM_DATA.u_no_sdslip = u_no_sdslip'; 
    SIM_DATA.v_no_sdslip = v_no_sdslip'; 
    SIM_DATA.w_no_sdslip = w_no_sdslip';
    SIM_DATA.V_no_sdslip = sqrt(u_no_sdslip.^2 + v_no_sdslip.^2 + w_no_sdslip.^2)';
    
    SIM_DATA.Lift_sdslip = -F_wind_z'; 
    SIM_DATA.Drag_sdslip = -F_wind_x';
    SIM_DATA.Normal_force = -F_ans_z'; 
    SIM_DATA.Axial_force = -F_ans_x';
    
    SIM_DATA.Fx_aero_b = Fx_aero_b'; 
    SIM_DATA.Fy_aero_b = Fy_aero_b'; 
    SIM_DATA.Fz_aero_b = Fz_aero_b';
    SIM_DATA.Fx_b = Fx_b'; 
    SIM_DATA.Fy_b = Fy_b'; 
    SIM_DATA.Fz_b = Fz_b';
    
    SIM_DATA.Mx_no_sdslip = Mx_no_sd'; 
    SIM_DATA.My_no_sdslip = My_no_sd'; 
    SIM_DATA.Mz_no_sdslip = Mz_no_sd';
    SIM_DATA.Mx_b = Mx_b'; 
    SIM_DATA.My_b = My_b'; 
    SIM_DATA.Mz_b = Mz_b';
    
    SIM_DATA.CL = C_L'; SIM_DATA.CD = C_D'; SIM_DATA.Cm = C_m';
    SIM_DATA.no_sdslip_angle = no_sdslip_angle'; 
    SIM_DATA.alpha_equivalent = alpha_eq';
    SIM_DATA.CoP_fraction = CoP_frac';
    
    % De-spun frame calculations for plots (requires roll angle phi)
    eul = quat2eul(x(:,7:10), 'ZYX');              
    phi = eul(:,3);
    c_despun = cos(phi); s_despun = sin(phi);
    
    SIM_DATA.Fx_nr =  Fx_b'; 
    SIM_DATA.Fy_nr =  Fy_b' .* c_despun + Fz_b' .* s_despun;
    SIM_DATA.Fz_nr = -Fy_b' .* s_despun + Fz_b' .* c_despun;
    SIM_DATA.Mx_nr =  Mx_b';
    SIM_DATA.My_nr =  My_b' .* c_despun + Mz_b' .* s_despun;
    SIM_DATA.Mz_nr = -My_b' .* s_despun + Mz_b' .* c_despun;

end