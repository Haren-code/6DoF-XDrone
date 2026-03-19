clear; clc; close all;

%% 1. Define Maneuver Parameters
V_body_x = 15;          % [m/s] Forward speed
Turn_Time = 10;         % [s] Duration for 180 deg turn
Turn_Angle = pi;        % [rad] 180 degrees
Omega_turn = Turn_Angle / Turn_Time; % [rad/s] Required Yaw Rate

dt = 0.01;              % Time step
T_end = Turn_Time;      % Simulation duration
time = 0:dt:T_end;

% Drone Params (Example)
m = 5.0;                % [kg] Mass
g = 9.81;               % [m/s^2] Gravity

%% 2. Calculate Equilibrium Bank Angle
% Centripetal Accel = V * Omega
a_lat = V_body_x * Omega_turn;

% Bank Angle (Phi) required to sustain this turn without slipping
phi_target = atan(a_lat / g);

fprintf('Maneuver Specs:\n');
fprintf('  Speed: %.1f m/s\n', V_body_x);
fprintf('  Yaw Rate: %.3f rad/s\n', Omega_turn);
fprintf('  Turn Radius: %.1f m\n', V_body_x / Omega_turn);
fprintf('  Required Bank Angle: %.1f deg\n', rad2deg(phi_target));

%% 3. Kinematic Simulation Loop
% State: [North, East, Altitude, Phi, Theta, Psi]
pos = [0; 0; 100];      % Start at [0,0,100]
euler = [phi_target; 0; 0]; % Start already banked [Roll, Pitch, Yaw]

path_log = zeros(3, length(time));

for i = 1:length(time)
    
    % Store Position
    path_log(:,i) = pos;
    
    % 1. Body Velocities (Constrained)
    u = V_body_x; 
    v = 0; % No side-slip (perfect coordinated turn)
    w = 0; % Constant altitude body frame approx
    
    % 2. Rotation Matrix (Body -> Inertial)
    phi = euler(1); theta = euler(2); psi = euler(3);
    
    R_z = [cos(psi) -sin(psi) 0; sin(psi) cos(psi) 0; 0 0 1];
    R_y = [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)];
    R_x = [1 0 0; 0 cos(phi) -sin(phi); 0 sin(phi) cos(phi)];
    R_ib = R_z * R_y * R_x;
    
    V_inertial = R_ib * [u; v; w];
    
    % 3. Update Rates (Kinematic Coordinated Turn)
    % Global yaw rate is constant Omega_turn.
    % We need to find body rates (p,q,r) that result in dPsi/dt = Omega_turn
    % Simple kinematic update for Euler angles:
    dot_euler = [0; 0; Omega_turn]; 
    
    % 4. Integration
    pos = pos + V_inertial * dt;
    euler = euler + dot_euler * dt;
    
end

%% 4. Visualization
figure('Color','w');
plot(path_log(2,:), path_log(1,:), 'b-', 'LineWidth', 2); grid on;
hold on;
plot(path_log(2,1), path_log(1,1), 'go', 'MarkerFaceColor','g'); % Start
plot(path_log(2,end), path_log(1,end), 'rs', 'MarkerFaceColor','r'); % End
axis equal;
xlabel('East [m]'); ylabel('North [m]');
title(['Coordinated 180^\circ Turn (R \approx ' num2str(V_body_x/Omega_turn, '%.1f') ' m)']);
subtitle(['Bank Angle: ' num2str(rad2deg(phi_target), '%.1f') '^\circ']);