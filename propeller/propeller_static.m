%% HQProp 3x3x3 Static Data Plotter
clear, clc

% 1. Load and Sort the Pre-Processed Dataset
filename = 'HQProp_3x3x3_Extracted.csv';
fprintf('Loading dataset: %s\n', filename);

% Use readmatrix to extract purely the numbers 
% The columns in the file are: 
% 1:RPM | 2:Thrust | 3:Watts | 4:C_T0 | 5:C_P0 | 6:C_Q0
data = readmatrix(filename, 'NumHeaderLines', 1);

% SORT THE ENTIRE DATASET BY RPM (Column 1)
data = sortrows(data, 1);

% Assign the now-sorted columns to simple arrays
rpm  = data(:, 1);
C_T0 = data(:, 4);
C_P0 = data(:, 5);

% 2. Calculate the Final Averages (Weighted by RPM)
% ensure  average RPMs > 10000
valid_idx = rpm > 10000;
rpm_valid = rpm(valid_idx);

% or using median or mean would work as well
% mean_CT0 = median(C_T0(valid_idx));
% mean_CP0 = median(C_P0(valid_idx));

mean_CT0 = sum(C_T0(valid_idx) .* rpm_valid) / sum(rpm_valid);
mean_CP0 = sum(C_P0(valid_idx) .* rpm_valid) / sum(rpm_valid);

fprintf('\n--- Extracted Static Coefficients ---\n');
fprintf('Mean C_T0 (Thrust): %.4f\n', mean_CT0);
fprintf('Mean C_P0 (Power Electric):  %.4f\n', mean_CP0);
fprintf('-------------------------------------\n');

% 3. Plot the Results
figure('Name', 'HQProp 3x3x3 Static Coefficients vs RPM', ...
       'Position', [100, 100, 900, 600]);

% Plot the data points (Now plotting cleanly from left to right!)
plot(rpm(valid_idx), C_T0(valid_idx), 'b-o', ...
    'LineWidth', 1.5, 'MarkerFaceColor', 'b', 'DisplayName', 'C_{T0} (Thrust Coeff)');
hold on;
plot(rpm(valid_idx), C_P0(valid_idx), 'r-o', ...
    'LineWidth', 1.5, 'MarkerFaceColor', 'r', 'DisplayName', 'C_{P0} (Power Coeff)');

% Plot the average dashed lines
yline(mean_CT0, 'b--', 'LineWidth', 1.5, 'HandleVisibility', 'off'); 
yline(mean_CP0, 'r--', 'LineWidth', 1.5, 'HandleVisibility', 'off');

% Formatting
title('HQProp 3x3x3 Static Coefficients vs RPM', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('RPM', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Coefficient Value', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 11);
grid on;