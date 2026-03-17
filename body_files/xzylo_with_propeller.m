function [properties,aerodynamics,initialization] = xzylo_with_propeller()
    %%--------------------- Parameters for inertia calculation ------------------------
    % Ring 1
    M1 = 6.23*1e-3;      % kg
    R_o1 = 97*1e-3 /2;   % m (Outer Radius)
    R_i1 = 96*1e-3 /2;   % m (Inner Radius)
    L1 = 54.5*1e-3;      % m (Length/Chord)
    X1 = 54.5*1e-3 /2;   % m (X-position of CoM)
    
    % Ring 2 
    M2 = 16.5*1e-3;      % kg
    R_o2 = 97*1e-3 /2;   % m
    R_i2 = 95*1e-3 /2;   % m
    L2 = 13*1e-3;        % m
    X2 = 13.0*1e-3 /2 ;  % m

    %% --------------- Execute the function for inertia calculation --------------------
    inertia_results = calculate_inertia(M1, R_o1, R_i1, L1, X1, M2, R_o2, R_i2, L2, X2);

    properties.Xzylo.CoG_pos = inertia_results.r_cg; % From the leading edge 0.0122
    properties.Xzylo.percentage_CoG = properties.Xzylo.CoG_pos./L1 ; %0.2236 original
    properties.Xzylo.I = inertia_results.I_total; %%% AROUND ITS OWN COG
    properties.Xzylo.mass = inertia_results.M_total;  % mass [kg]
    properties.Xzylo.b = R_o1;   % span = diameter
    properties.Xzylo.c = L1;     % chord = distance from LE to TE
    properties.Xzylo.S = properties.Xzylo.b*properties.Xzylo.c; % reference area
    properties.Xzylo.aspect_ratio = properties.Xzylo.b/properties.Xzylo.c; % span/chord
    properties.Xzylo.f_coeff = 5e-3; % Friction Torque Coefficient

    %% ------------------------- PROPELLER CONFIGURATION ---------------------------

    % Negative RPM = Counter-Rotating (Opposite to X-Zylo spin)
    Prop_rpm   = -25000;  % 10 - 20k rpm 
    properties.Prop.omega = 2*pi/60 * Prop_rpm;
    properties.Prop.CoG_pos = [0.0 ; 0; 0] + properties.Xzylo.CoG_pos;    % Position relative to CoG [m] (e.g. [0.05;0;0] for forward)
    properties.Prop.mass  = 0.3 *1e-3;          % Mass of motor + prop + mount [kg] (Estimate 0.3g)
    properties.Prop.radius = 5.5/2 *1e-3; 

    % 40mm diameter K_t = 6.8*10-9 N/(rad/s)^2
    % k_q = 0.02 * torque
    % 67mm dimatere K_t = 2.66*e-7
    properties.Prop.k_T   = 4e-8;                % Thrust coeff [N / Rad/s^2] (Estimate)
    properties.Prop.k_V   = 0;
    properties.Prop.k_Q   = 0.5*properties.Prop.radius*properties.Prop.k_T; % 2% of the thrust one OR 0.1*diameter*k_t
     
    % Motor
    properties.Motor.mass = 3.54*1e-3;      % kg
    properties.Motor.radius = 7/2*1e-3;   % m (Radius)
    properties.Motor.L = 20e-3;      % m (Length/Chord)
    properties.Motor.CoG_pos = properties.Prop.CoG_pos + [properties.Motor.L/2; 0; 0]; %position of COG of the motor

    Motor.I_xx = 0.5*properties.Motor.mass*properties.Motor.radius^2;
    Motor.I_yyzz = (0.5*Motor.I_xx + 1/12*properties.Motor.mass*properties.Motor.L^2);
    

    %% ------------------------- Complete System Mass and CoG -------------------------
    properties.mass_total = properties.Motor.mass + properties.Prop.mass + properties.Xzylo.mass;
    properties.CoG_pos_total = (properties.Xzylo.mass*properties.Xzylo.CoG_pos + properties.Motor.mass*properties.Motor.CoG_pos + properties.Prop.mass*properties.Prop.CoG_pos)/properties.mass_total;
    properties.percentage_CoG_total = properties.CoG_pos_total./L1;

    % X-zylo MoI
    properties.Xzylo.I(2,2) = properties.Xzylo.I(2,2) + properties.Xzylo.mass*(properties.Xzylo.CoG_pos(1)-properties.CoG_pos_total(1))^2;
    properties.Xzylo.I(3,3) = properties.Xzylo.I(2,2);

    % Propeller MoI
    prop.I_xx = 0.5*properties.Prop.mass*properties.Prop.radius^2;
    properties.Prop.I_yyzz = 0.5*prop.I_xx + properties.Prop.mass*(properties.Prop.CoG_pos(1)-properties.CoG_pos_total(1))^2;
    properties.Prop.I = diag([prop.I_xx , properties.Prop.I_yyzz, properties.Prop.I_yyzz]); 

    % Motor MoI
    Motor.I_xx = 0.5*properties.Motor.mass*properties.Motor.radius^2;
    Motor.I_yyzz_0 = (0.5*Motor.I_xx + 1/12*properties.Motor.mass*properties.Motor.L^2);
    Motor.I_yyzz = Motor.I_yyzz_0 + properties.Motor.mass*(properties.Motor.CoG_pos(1)-properties.CoG_pos_total(1))^2;
    properties.Motor.I = diag([Motor.I_xx , Motor.I_yyzz, Motor.I_yyzz]); 


    %% Aerodynamics Group
    % ------------------------- Environmental Properties ---------------------------
    aerodynamics.rho = 1.225;
    aerodynamics.V_wind_i = [0; 0; 0]; % wind velocity in inertial frame [m/s]
    aerodynamics.g = 9.81;
    
    %----------------------- Aerodynamic coefficient functions --------------------------
    aerodynamics_xzylo      % run aerodynamic file for interpolation of data
    aerodynamics.Xzylo.C_L =  @(angle) interp1(alpha_rad_CL, C_L_sym_smooth, angle, 'pchip', 'extrap');            
    aerodynamics.Xzylo.C_D =  @(angle) interp1(alpha_rad_CD, C_D_sym_smooth, angle, 'pchip', 'extrap');
    aerodynamics.Xzylo.CoP_frac = @(angle) interp1(alpha_rad_COP, COP_sym_smooth, angle, 'pchip', 'extrap')./100;
    aerodynamics.Xzylo.C_m =  @(angle) (properties.percentage_CoG_total(1) - aerodynamics.Xzylo.CoP_frac(angle)) * (aerodynamics.Xzylo.C_L(angle) * cos(angle) + aerodynamics.Xzylo.C_D(angle) * sin(angle));
    aerodynamics.Xzylo.C_Y =  @(angle) 0;

    % Optional external force/moment (written in INERTIAL frame)
    aerodynamics.Fext = @(t) (t >= 0 && t <= 30) * [0; 0; 0];
    aerodynamics.Mext = @(t) (t >= 0 && t <= 30) * [0; 0; 0]; 

    % ------------------------- Initialization states ------------------------------------
    initialization.launch_angle = 10; % launch angle in degrees
    initialization.V_mag = 20; % Magnitude of the launch velocity
    initialization.Omega_mag = 50; % Rotational speed at launch [RPS]

    initialization.v0 = [initialization.V_mag*cosd(0),  0  , -initialization.V_mag*sind(0)]; % initial velocity vector (u,v,w)_0   [m/s]
    initialization.omega0 = [2*pi*initialization.Omega_mag, 0, 0]; % initial rotational velocity vector (p,q,r)_0  [rad/s]
    initialization.euler0 = [0*pi/180, initialization.launch_angle*pi/180, 0*pi/180]; %[yaw pitch roll];   % [psi theta phi]
    initialization.quat0  = eul2quat(initialization.euler0);   % returns [w x y z]
    initialization.pos0 = [0 0 -30]; % initial position in inertial frame [m]
    initialization.tf = 30; %maximum simulation time
   
end