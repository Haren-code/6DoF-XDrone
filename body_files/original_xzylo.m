function [properties,aerodynamics,initialization] = original_xzylo()

    %% --------------- Properties of the craft - X-Zylo in this case -------------------
    properties.Xzylo.CoG_pos = 12.234e-3; % From the leading edge 
    properties.Xzylo.b = 96.5e-3;         % span = diameter
    properties.Xzylo.c = 54.5e-3;         % chord = distance from LE to TE
    properties.Xzylo.percentage_CoG =  properties.Xzylo.CoG_pos/properties.Xzylo.c; %0.22448
    properties.Xzylo.mass = 22.69e-3;  % mass [kg]
    properties.Xzylo.I = 1e-9 * [52309.988 0 0
                                 0 29891.828 0
                                 0 0 29891.828]; % Around the CoG
    properties.Xzylo.S = properties.Xzylo.b*properties.Xzylo.c; % reference area
    properties.Xzylo.aspect_ratio = properties.Xzylo.b/properties.Xzylo.c; % span/chord


    %% ------------------------- Complete System Mass and CoG -------------------------
    properties.mass_total = properties.Xzylo.mass;
    properties.CoG_pos_total = properties.Xzylo.CoG_pos;
    properties.percentage_CoG_total = properties.Xzylo.percentage_CoG;
    properties.Area = properties.Xzylo.S;
    properties.chord = properties.Xzylo.c;
    properties.span = properties.Xzylo.b;
    properties.f_coeff = 7e-4; % Friction Torque Coefficient
    
    % ------------------------- Environmental Properties ---------------------------
    properties.rho = 1.225;
    properties.V_wind_i = [0; 0; 0]; % wind velocity in inertial frame [m/s]
    properties.g = 9.81;
       
    %% Aerodynamics Group

    %----------------------- Aerodynamic coefficient functions --------------------------
    aerodynamics_xzylo      % run aerodynamic file for interpolation of data    
    aerodynamics.Xzylo.C_L = @(angle) C_L_interp(angle);
    aerodynamics.Xzylo.C_D = @(angle) C_D_interp(angle);
    aerodynamics.Xzylo.CoP_frac = @(angle) CoP_interp(angle);
    % aerodynamics.Xzylo.C_m =  @(angle) (x_CoG - x_CoP(angle)) * (C_L*sin(angle) + C_D*cos(angle))
    % aerodynamics.Xzylo.C_Y =  @(angle) 0;

    % Trim conditions 
    [properties.alpha_trim, properties.V_trim, properties.Thrust_req] = calculate_trim(properties.percentage_CoG_total, properties.mass_total*properties.g, properties.rho,...
                                                                                    properties.Area, aerodynamics.Xzylo.C_L, aerodynamics.Xzylo.C_D, aerodynamics.Xzylo.CoP_frac);

    % Optional external force/moment (written in BODY frame)
    aerodynamics.Fext = @(t) (t >= 0 && t <= 30) * [0; 0; 0];
    aerodynamics.Mext = @(t) (t >= 0 && t <= 30) * [0; 0; 0]; 

    % ------------------------- Initialization states ------------------------------------
    initialization.launch_angle = 12; % launch angle in degrees
    initialization.V_mag = 20; % Magnitude of the launch velocity
    initialization.Omega_mag = 40; % Rotational speed at launch [RPS]
    initialization.AoA = 0;    

    initialization.v0 = [initialization.V_mag*cosd(initialization.AoA),  0  , initialization.V_mag*sind(initialization.AoA)]; % initial velocity vector (u,v,w)_0   [m/s]
    initialization.omega0 = [2*pi*initialization.Omega_mag, 0, 0]; % initial rotational velocity vector (p,q,r)_0  [rad/s]
    initialization.euler0 = [deg2rad(0), deg2rad(initialization.launch_angle), deg2rad(0)]; %[yaw pitch roll];   % [psi theta phi]
    initialization.quat0  = eul2quat(initialization.euler0);   % returns [w x y z]
    initialization.pos0 = [0 0 -1.5]; % initial position in inertial frame [m]
    initialization.tf = 30; %maximum simulation time
   
end