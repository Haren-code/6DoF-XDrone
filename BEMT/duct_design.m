clear; clc; close all;
%% ========================================================================
%%                    PART 1: DUCT GEOMETRY GENERATION
%% ========================================================================
%% ============ USER INPUTS ============
csv_filename  = 'xf-n0012-il-200000.csv'; 
geom_filename = 'NACA0012.dat';           

% -- Duct Dimensions --
R_duct_mid = 56.7e-3; % [m] Radius of the duct ring (centerline)
duct_chord = 0.1;     % [m] Length of the duct (Chord)
N_duct     = 60;      % Azimuthal resolution (strip theory)
rho        = 1.225;

%% -------- DUCT GEOMETRY LOADING (FROM NACA0012.dat) --------
fprintf('--- Loading Duct Geometry from %s ---\n', geom_filename);
if isfile(geom_filename)
    % Logic matched to your original code
    raw_data = readmatrix(geom_filename, 'NumHeaderLines', 2);
    
    % Find the split between Upper and Lower surfaces
    split_idx = find(diff(raw_data(:,1)) < -0.5, 1);
    
    if isempty(split_idx)
        error('Could not identify Upper/Lower split in airfoil file.');
    else
        Upper_X = raw_data(1:split_idx, 1); Upper_Y = raw_data(1:split_idx, 2);
        Lower_X = raw_data(split_idx+1:end, 1); Lower_Y = raw_data(split_idx+1:end, 2);
        
        % Stitch them together (Upper -> TE -> Lower)
        X_segs = [Upper_X; 1.0; flipud(Lower_X)];
        Y_segs = [Upper_Y; 0.0; flipud(Lower_Y)];
        X_norm = X_segs'; 
        Y_norm = Y_segs';
    end
    fprintf('Loaded %s successfully.\n', geom_filename);
else
    error('Airfoil file %s not found!', geom_filename);
end

% Scale Airfoil (Chord along X, Thickness along Y)
X_af = X_norm * duct_chord;
Y_af = Y_norm * duct_chord;

%% ============== FIGURE 1: 3D DUCT GEOMETRY (RESTORED) ==================
fig = figure('Color','w','Name','Duct Geometry');
axes; hold on; axis equal; grid on; view(30, 20);
xlabel('X (Axial) [m]'); ylabel('Y [m]'); zlabel('Z [m]');
title(['Duct Geometry (Ring Wing) - ' geom_filename]);

% Construct 3D Surface
theta = linspace(0, 2*pi, N_duct);
X_3D = []; Y_3D = []; Z_3D = [];
for k = 1:length(theta)
    th = theta(k);
    % Rotate Airfoil Y (Thickness) into Radial Direction
    % X remains X (Axial)
    % Radial Offset = R_duct_mid
    
    R_local = R_duct_mid + Y_af; % Radial position of surface points
    
    x_ring = X_af;              % Axial coordinate
    y_ring = R_local * cos(th); % Y coordinate
    z_ring = R_local * sin(th); % Z coordinate
    
    X_3D = [X_3D; x_ring];
    Y_3D = [Y_3D; y_ring];
    Z_3D = [Z_3D; z_ring];
end
surf(X_3D, Y_3D, Z_3D, 'FaceColor',[0.6 0.6 0.8], 'EdgeColor','none', 'FaceAlpha', 0.8);
camlight; lighting gouraud;
plot3([-0.05 0.15], [0 0], [0 0], 'r-', 'LineWidth', 2); % Axis Marker

%% ========================================================================
%%            PART 2: DUCT AERODYNAMICS (STRIP THEORY)
%% ========================================================================
fprintf('\n--- Calculating Duct Forces ---\n');

% 1. Load Polar Data (Robust)
if isfile(csv_filename)
    opts = detectImportOptions(csv_filename);
    opts.VariableNamesLine = 11; opts.DataLines = [12 Inf];
    PolarData = readtable(csv_filename, opts);
    
    vars = PolarData.Properties.VariableNames;
    idx_a = find(strcmpi(vars, 'Alpha'), 1);
    idx_cl = find(strcmpi(vars, 'Cl'), 1);
    idx_cd = find(strcmpi(vars, 'Cd'), 1);
    
    Alpha_dat = PolarData{:, idx_a};
    Cl_dat    = PolarData{:, idx_cl};
    Cd_dat    = PolarData{:, idx_cd};
    [Alpha_dat, uniqueIdx] = unique(Alpha_dat);
    Cl_dat = Cl_dat(uniqueIdx); Cd_dat = Cd_dat(uniqueIdx);
    fprintf('Loaded Polar Data: %s\n', csv_filename);
else
    error('CSV file %s not found!', csv_filename);
end

% 2. Simulation Grid
alpha_vec = linspace(-45, 45, 91); % Restrict to -45 to +45
v_vec     = linspace(0, 30, 10);
[ALPHA_GRID, V_GRID] = meshgrid(alpha_vec, v_vec);

Lift_Map = zeros(size(ALPHA_GRID));
Drag_Map = zeros(size(ALPHA_GRID));

% Ring Properties
circumference = 2 * pi * R_duct_mid;
strip_width   = circumference / N_duct;
strip_area    = strip_width * duct_chord;

% =========================================================================
%                  PCHIP SPLINE SETUP
% =========================================================================
% We add points at +/- 45 degrees.
% PCHIP prevents the curve from shooting up (overshooting).
Alpha_Aug = [ -45;      Alpha_dat;   45    ];
Cl_Aug    = [ 0;        Cl_dat;      0     ];
Cd_Aug    = [ 1.2;      Cd_dat;      1.2   ];

F_cl_final = griddedInterpolant(Alpha_Aug, Cl_Aug, 'pchip', 'nearest');
F_cd_final = griddedInterpolant(Alpha_Aug, Cd_Aug, 'pchip', 'nearest');
% =========================================================================

for i = 1:size(ALPHA_GRID, 1)     % Velocity Loop
    for j = 1:size(ALPHA_GRID, 2) % Alpha Loop
        
        V_inf = V_GRID(i,j);
        AoA_glob = deg2rad(ALPHA_GRID(i,j));
        
        V_x_glob = V_inf * cos(AoA_glob); 
        V_z_glob = V_inf * sin(AoA_glob); 
        
        F_x_total = 0; 
        F_z_total = 0; 
        
        % --- STRIP THEORY INTEGRATION ---
        for th = linspace(0, 2*pi, N_duct)
            
            U_c = V_x_glob; 
            U_n = V_z_glob * sin(th); 
            
            % Local Angle of Attack
            if abs(U_c) < 1e-3
                 local_alpha_deg = 90 * sign(U_n);
                 V_local = abs(U_n);
            else
                 local_alpha_deg = rad2deg(atan2(U_n, U_c));
                 V_local = sqrt(U_c^2 + U_n^2);
            end
            
            % --- LOOKUP (Spline with Clamping to +/- 45) ---
            q_alpha = max(min(local_alpha_deg, 45), -45);
            cl = F_cl_final(q_alpha);
            cd = F_cd_final(q_alpha);
            % -----------------------------------------------
            
            % Forces in Strip Frame
            q = 0.5 * rho * V_local^2 * strip_area;
            L_strip = q * cl;
            D_strip = q * cd;
            
            gamma = atan2(U_n, U_c);
            
            F_x_strip = L_strip * (-sin(gamma)) + D_strip * (-cos(gamma)); 
            F_r_strip = L_strip * (cos(gamma))  + D_strip * (-sin(gamma));
            
            F_z_strip = F_r_strip * sin(th);
            
            F_x_total = F_x_total + F_x_strip;
            F_z_total = F_z_total + F_z_strip;
        end
        
        Drag_Map(i,j) = -(F_x_total * cos(AoA_glob) + F_z_total * sin(AoA_glob));
        Lift_Map(i,j) =  (-F_x_total * sin(AoA_glob) + F_z_total * cos(AoA_glob));
        
    end
end
fprintf('Done.\n');

%% ============== FIGURE 2: FORCE MAPS ==================
figure('Color','w','Name','Duct Aerodynamics','WindowStyle','docked');
subplot(1,2,1); 
contourf(ALPHA_GRID, V_GRID, Lift_Map, 20); colorbar;
xlabel('Drone AoA [deg]'); ylabel('Velocity [m/s]');
title('Duct Net Lift [N]'); grid on;

subplot(1,2,2); 
contourf(ALPHA_GRID, V_GRID, Drag_Map, 20); colorbar;
xlabel('Drone AoA [deg]'); ylabel('Velocity [m/s]');
title('Duct Net Drag [N]'); grid on;

%% ============== FIGURE 3: 1D CUT ==================
figure('Color','w','Name','Duct 1D Polars');
mid_idx = floor(size(V_GRID,1)/2); 
V_plot = V_GRID(mid_idx, 1);

plot(alpha_vec, Lift_Map(mid_idx, :), 'b-o', 'LineWidth', 2, 'DisplayName', 'Lift');
hold on;
plot(alpha_vec, Drag_Map(mid_idx, :), 'r-s', 'LineWidth', 2, 'DisplayName', 'Drag');
grid on; legend;
xlabel('Drone Angle of Attack [deg]'); ylabel('Force [N]');
title(['Duct Forces @ ' num2str(V_plot) ' m/s']);
xlim([-45 45]);

%% ============== FIGURE 4: DEBUG SPLINE ==================
figure('Color','w','Name','DEBUG: PCHIP Check');

a_check = linspace(-45, 45, 200);
cl_check = F_cl_final(a_check);
cd_check = F_cd_final(a_check);

subplot(2,1,1); hold on; grid on;
plot(Alpha_dat, Cl_dat, 'ko', 'MarkerFaceColor','k', 'DisplayName','Raw Data');
plot(a_check, cl_check, 'b-', 'LineWidth', 2, 'DisplayName','PCHIP Fit');
ylabel('Cl'); title('Lift Coeff (Shape Preserving)'); 
legend('Location','best'); xlim([-45 45]);

subplot(2,1,2); hold on; grid on;
plot(Alpha_dat, Cd_dat, 'ko', 'MarkerFaceColor','k', 'DisplayName','Raw Data');
plot(a_check, cd_check, 'r-', 'LineWidth', 2, 'DisplayName','PCHIP Fit');
ylabel('Cd'); xlabel('Alpha [deg]'); title('Drag Coeff'); 
xlim([-45 45]);