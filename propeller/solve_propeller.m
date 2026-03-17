function [J, RPM, Torque, C_P] = solve_propeller(V, Thrust_req)
% SOLVE_PROPELLER Calculates required RPM, Advance Ratio, and Mechanical Power
%
% INPUTS:
%   V          - Flight velocity (m/s)
%   Thrust_req - Required thrust per motor (Newtons)
%
% OUTPUTS:
%   J          - Advance Ratio
%   RPM        - Required motor speed
%   P_mech     - Mechanical shaft power required (Watts)

    %% 1. Propeller & Environment Constants
    D_in = 3.0;             % Diameter in inches
    pitch_in = 3.0;         % Pitch in inches
    D_m = D_in * 0.0254;    % Convert diameter to meters
    
    rho = 1.225;            % Air density (kg/m^3)
    
    % Aerodynamic properties from your static bench data
    C_T0 = 0.205;           
    C_P0 = C_T0 / 1.5;      
    unload_factor = 0.7;    
    
    PD_ratio = pitch_in / D_in;
    J_max = PD_ratio * 0.955;

    %% 2. Handle Zero-Velocity (Hover) Condition
    % If V=0, J=0. We can solve for n algebraically without the numerical solver.
    if V == 0
        n_rps = sqrt(Thrust_req / (C_T0 * rho * D_m^4));
        
    %% 3. Handle Forward Flight (Numerical Solver)
    else
        % Define the nonlinear thrust equation to solve for n (revolutions per sec)
        % Equation format: Thrust_calculated(n) - Thrust_required = 0
        thrust_eq = @(n) (C_T0 .* (1 - ( (V ./ (n .* D_m .* J_max)) ).^1.5)) .* rho .* n.^2 .* D_m.^4 - Thrust_req;
        
        % Generate a smart initial guess for the solver based on hover RPM
        n_guess = sqrt(Thrust_req / (C_T0 * rho * D_m^4));
        
        % Ensure the guess is fast enough to prevent negative J terms (windmilling)
        n_min = V / (J_max * D_m);
        if n_guess <= n_min
            n_guess = n_min * 1.05; 
        end
        
        % Run the numerical solver to find the exact 'n'
        options = optimset('Display', 'off'); % Run silently
        n_rps = fzero(thrust_eq, n_guess, options);
    end

    %% 4. Calculate Final Outputs
    RPM = n_rps * 60;
    
    if n_rps > 0
        J = V / (n_rps * D_m);
    else
        J = 0;
    end
    
    % Calculate the Aerodynamic Power Coefficient at this J
    C_P = C_P0 * (1 - (unload_factor * (J / J_max)^1.33));

    % Calculate the Torque Coefficient (C_Q)
    C_Q = C_P / (2 * pi);

    % Calculate final Mechanical Power in Watts
    Power_mech = C_P * rho * n_rps^3 * D_m^5;
    
    % Calculate final Aerodynamic Torque (assuming standard SI units, this will be in Nm)
    Torque = C_Q * rho * n_rps^2 * D_m^5;

end