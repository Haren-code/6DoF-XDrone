%% Propeller Aerodynamic Coefficient Predictor
% Models C_T, C_P, C_Q, and Efficiency vs Advance Ratio (J)
% Tuned for an aggressive 3x3x3 (Square) Multicopter Propeller

clear; close all; clc;

%% 1. Propeller Geometry and Static Bench Inputs
D_in = 3.0;       % Diameter in inches
pitch_in = 3.0;   % Pitch in inches
blades = 3;       % Number of blades

% --- STATIC BENCH DATA (J = 0) ---
% Replace these with your actual calculated test bench values!
C_T0 = 0.205;         % Estimated static thrust coefficient 
C_P0 = C_T0/1.5;      % Estimated static power coefficient 

%% 2. Flight Dynamic Model Parameters
% Calculate Pitch-to-Diameter Ratio
PD_ratio = pitch_in / D_in;

% Estimate the Zero-Thrust Advance Ratio (J_max)
% Airfoils typically stop lifting slightly past their geometric pitch 
J_max = PD_ratio * 0.955; 

% Unload Factor: How much does the power drop at top speed?
% A square prop (P/D = 1.0) unloads heavily. We assume a 70% drop.
unload_factor = 0.7; 

%% 3. Generate the Advance Ratio Array (J)
% Sweep J from 0 to slightly past J_max to see the zero-crossing
J = linspace(0, J_max * 1.25, 100);

%% 4. Calculate Dynamic Coefficients
% Thrust Coefficient (C_T): Modeled as power drop to zero at J_max
C_T = C_T0 .* (1 - (J ./ J_max).^1.5);

% Power Coefficient (C_P): Drops by the unload factor at J_max
C_P = C_P0 .* (1 - (unload_factor .* (J ./ J_max).^1.33));

% Torque Coefficient (C_Q): Directly derived from Power Coefficient
C_Q = C_P ./ (2 * pi);

% Efficiency (eta): Calculate and convert to percentage
eta = J .* (C_T ./ C_P);
eta_percent = eta .* 100;

% Physical constraint: Efficiency cannot be negative (if J > J_max, prop is braking)
eta_percent(C_T < 0) = 0; 

%% 5. Plotting the Results
figure('Name', '3x3x3 Propeller Performance Curves', 'Position', [100, 100, 900, 600]);

% Plot 1: Thrust Coefficient (C_T)
subplot(2,2,1);
plot(J, C_T, 'b-', 'LineWidth', 2);
yline(0, 'k--'); % Zero line
title('Thrust Coefficient (C_T)');
xlabel('Advance Ratio (J)');
ylabel('C_T');
grid on;
xlim([0, J_max]);      % Restricts the X-axis from 0 to J_max
ylim([0, inf]);        % Restricts the Y-axis to start at 0, 'inf' lets MATLAB auto-scale the top

% Plot 2: Power & Torque Coefficients (C_P & C_Q)
% Plotting together since they share the exact same curve shape, just scaled
subplot(2,2,2);
yyaxis left
plot(J, C_P, 'r-', 'LineWidth', 2);
ylabel('Power Coefficient (C_P)', 'Color', 'r');
set(gca, 'YColor', 'r');

yyaxis right
plot(J, C_Q, 'm--', 'LineWidth', 1.5);
ylabel('Torque Coefficient (C_Q)', 'Color', 'm');
set(gca, 'YColor', 'm');

title('Power & Torque Coefficients');
xlabel('Advance Ratio (J)');
grid on;

% Plot 3: Propeller Efficiency (eta)
subplot(2,2,[3,4]); % Spans the bottom row
plot(J, eta_percent, 'g-', 'LineWidth', 2.5);
title('Propeller Aerodynamic Efficiency (\eta)');
xlabel('Advance Ratio (J)');
ylabel('Efficiency (%)');
ylim([0, max(eta_percent) * 1.2]); % Scale y-axis nicely
grid on;

% Add a marker at peak efficiency
[max_eta, max_idx] = max(eta_percent);
hold on;
plot(J(max_idx), max_eta, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
text(J(max_idx) + 0.05, max_eta, sprintf('Peak: %.1f%% at J=%.2f', max_eta, J(max_idx)), ...
    'FontWeight', 'bold', 'FontSize', 10);