function [properties,aerodynamics,initialization] = xzylo_with_propeller()
  
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