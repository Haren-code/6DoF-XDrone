%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Aerodynamic Analysis - Function only for plotting and analysis
%
% Outputs:
%   - CL, CD, CoP (raw + smooth + interpolated)
%   - Pitching moment about user-defined CG
%   - Static stability derivative
%   - Trim angle (Cm = 0)
%   - Zero stability point (dCm/dα = 0)
%   - Level flight Trim Speed and Required Thrust
%
% Convention:
%   Positive Cm = Nose-up
%   Stability requires dCm/dα < 0
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars; close all; clc;

%% ========================================================================
% USER INPUT
% ========================================================================
rho = 1.225; % Air density at sea level [kg/m^3]
g = 9.81;    % Gravity [m/s^2]

% --- Vehicle Selection ---
fprintf('\n========== VEHICLE SELECTION ==========\n');
fprintf('  1: X-Zylo  (x_CoG=22.85%% c, m=0.02273 kg, S=0.0052865 m^2)\n');
fprintf('  2: X-Drone (x_CoG=28.40%% c, m=0.1413 kg,  S=0.014375 m^2)\n');
fprintf('  3: Custom\n');
vehicle_choice = input('Select option [1/2/3]: ');

switch vehicle_choice
    case 1
        x_CoG   = 0.22448;   % midpoint of 22.4-23.3%
        m_drone = 0.02269;
        S_ref   = 0.00525925;
        fprintf('>> X-Zylo selected.\n');
    case 2
        x_CoG   = 0.284;
        m_drone = 0.1413;
        S_ref   = 0.014375;
        fprintf('>> X-Drone selected.\n');
    case 3
        fprintf('\nEnter CG location as fraction of chord (e.g. 0.25)\n');
        x_CoG = input('x_CoG = ');
        if ~isscalar(x_CoG) || x_CoG < 0 || x_CoG > 1
            error('x_CoG must be between 0 and 1.');
        end
        fprintf('Enter drone mass in kg (e.g. 1.5)\n');
        m_drone = input('mass (kg) = ');
        if ~isscalar(m_drone) || m_drone <= 0
            error('Mass must be a positive scalar.');
        end
        fprintf('Enter reference area (S) in m^2 (e.g. 0.12)\n');
        S_ref = input('S (m^2) = ');
        if ~isscalar(S_ref) || S_ref <= 0
            error('S_ref must be a positive scalar.');
        end
    otherwise
        error('Invalid vehicle selection. Choose 1, 2, or 3.');
end

% --- Aerodynamics File Selection ---
fprintf('\n======= AERODYNAMICS FILE =============\n');
fprintf('  1: aerodynamics_xzylo (default)\n');
fprintf('  2: Custom file\n');
aero_choice = input('Select option [1/2]: ');

switch aero_choice
    case 1
        aero_file = 'aerodynamics_xzylo';
    case 2
        % Look in the same folder as this script
        script_dir = fileparts(mfilename('fullpath'));
        m_files = dir(fullfile(script_dir, '*.m'));
        
        % Filter out this script itself
        this_script = [mfilename, '.m'];
        m_files = m_files(~strcmp({m_files.name}, this_script));
        
        if isempty(m_files)
            error('No other .m files found in %s', script_dir);
        end
        
        fprintf('\nAvailable .m files in script directory:\n');
        for k = 1:numel(m_files)
            fprintf('  %d: %s\n', k, m_files(k).name);
        end
        file_idx = input(sprintf('Select file [1-%d]: ', numel(m_files)));
        if file_idx < 1 || file_idx > numel(m_files)
            error('Invalid selection.');
        end
        [~, aero_file, ~] = fileparts(m_files(file_idx).name); % strip .m extension
    otherwise
        error('Invalid aerodynamics selection. Choose 1 or 2.');
end

fprintf('>> Loading: %s.m\n', aero_file);

%% ========================================================================
% LOAD AERODYNAMIC DATA
% ========================================================================
run(aero_file);

%% ========================================================================
% INTERPOLATION
% ========================================================================
alpha_range = deg2rad(linspace(0, 90, 1000));

C_L_fit = C_L_interp(alpha_range);
C_D_fit = C_D_interp(alpha_range);
CoP_fit = CoP_interp(alpha_range);

%% ========================================================================
% FORCE AND MOMENT COMPUTATION
% ========================================================================

% Normal force coefficient
C_N = C_L_fit .* cos(alpha_range) + C_D_fit .* sin(alpha_range);

% Correct sign convention:
% Positive Cm = Nose-up
C_m = (x_CoG - CoP_fit./100) .* C_N;

%% ========================================================================
% STABILITY DERIVATIVE
% ========================================================================

dCm_dAlpha = gradient(C_m, alpha_range);

%% ========================================================================
% TRIM SOLUTION (Cm = 0)
% ========================================================================

[alpha_trim, V_trim, Thrust_req, C_L_trim, C_D_trim] = calculate_trim(x_CoG, m_drone*g, rho, S_ref, C_L_interp, C_D_interp, CoP_interp);


%% ========================================================================
% ZERO GRADIENT POINT (Neutral Stability)
% ========================================================================

dCm_fun = @(a) interp1(alpha_range, dCm_dAlpha, a, 'pchip');

alpha_interval = [deg2rad(0.1), deg2rad(45)]; % degrees, choose interval around expected trim

alpha_neutral = NaN;
try
    alpha_neutral = fzero(dCm_fun, alpha_interval);
catch
    warning('No zero-gradient found in interval [%g, %g].', alpha_interval(1), alpha_interval(2));
end


%% ========================================================================
% PRINT RESULTS
% ========================================================================

fprintf('\n================ RESULTS ================\n');

if ~isnan(alpha_trim)
    fprintf('Trim angle: %.3f deg\n', rad2deg(alpha_trim));
    if C_L_trim > 0
        fprintf('Trim C_L: %.4f\n', C_L_trim);
        fprintf('Trim C_D: %.4f\n', C_D_trim);
        fprintf('Req. Flight Speed: %.2f m/s (%.1f km/h)\n', V_trim, V_trim * 3.6);
        fprintf('Req. Thrust : %.2f N\n', Thrust_req);
    else
        fprintf('Req. Flight Speed: N/A (C_L at trim is <= 0)\n');
    end
else
    fprintf('Trim angle: Not found\n');
end

if ~isnan(alpha_neutral)
    fprintf('Zero stability point: %.3f deg\n', rad2deg(alpha_neutral));
else
    fprintf('Zero stability point: Not found\n');
end

[~, idx0] = min(abs(alpha_range));
fprintf('dCm/dα at α≈0: %.6f 1/rad\n', dCm_dAlpha(idx0));

if dCm_dAlpha(idx0) < 0
    fprintf('Configuration is statically STABLE near α=0\n');
else
    fprintf('Configuration is statically UNSTABLE near α=0\n');
end

fprintf('=========================================\n\n');

%% ========================================================================
% PLOTS
% ========================================================================

% --- Lift ---
figure;
plot(rad2deg(alpha_rad_CL), C_L,'o','DisplayName','Raw'); hold on;
plot(rad2deg(alpha_rad_CL), C_L_smooth,'g--','LineWidth',1.5,'DisplayName','Smoothed');
plot(rad2deg(alpha_range), C_L_fit,'r-','LineWidth',1.5,'DisplayName','Interpolated');
xlabel('\alpha [deg]'); ylabel('C_L');
title('Lift Curve');
legend('Location','best'); grid on;

% --- Drag ---
figure;
plot(rad2deg(alpha_rad_CD), C_D,'o','DisplayName','Raw'); hold on;
plot(rad2deg(alpha_rad_CD), C_D_smooth,'g--','LineWidth',1.5,'DisplayName','Smoothed');
plot(rad2deg(alpha_range), C_D_fit,'r-','LineWidth',1.5,'DisplayName','Interpolated');
xlabel('\alpha [deg]'); ylabel('C_D');
title('Drag Curve');
legend('Location','best'); grid on;

% --- CoP ---
figure;
plot(rad2deg(alpha_rad_CoP), CoP,'o','DisplayName','Raw'); hold on;
plot(rad2deg(alpha_rad_CoP), CoP_smooth,'g--','LineWidth',1.5,'DisplayName','Smoothed');
plot(rad2deg(alpha_range), CoP_fit,'r-','LineWidth',1.5,'DisplayName','Interpolated');
xlabel('\alpha [deg]'); ylabel('CoP (% chord)');
title('Center of Pressure');
legend('Location','best'); grid on;

% --- Pitching Moment ---
figure;
plot(rad2deg(alpha_range), C_m,'b-','LineWidth',1.5); hold on;
yline(0,'k--');

if ~isnan(alpha_trim)
    plot(rad2deg(alpha_trim),0,'ro','MarkerSize',8,'LineWidth',2);
end

xlabel('\alpha [deg]');
ylabel('C_m');
title(sprintf('Pitching Moment about x_{CG}=%.2f c',x_CoG));
grid on;

% --- Stability Derivative ---
figure;
plot(rad2deg(alpha_range), dCm_dAlpha,'m-','LineWidth',1.5); hold on;
yline(0,'k--');

if ~isnan(alpha_neutral)
    plot(rad2deg(alpha_neutral),0,'ro','MarkerSize',8,'LineWidth',2);
end

xlabel('\alpha [deg]');
ylabel('dC_m/d\alpha [1/rad]');
title('Static Stability Derivative');
grid on;