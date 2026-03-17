function [alpha_trim, V_trim, Thrust_req, CL_t, CD_t] = calculate_trim(x_CoG, W, rho, S_ref, C_L_interp, C_D_interp, COP_interp)
    
    % 1. Define Pitching Moment Function (Cm)
    % Positive Cm = Nose-up
    Cm_fun = @(a) (x_CoG - (COP_interp(a)./100)) .* ...
                 (C_L_interp(a).*cos(a) + C_D_interp(a).*sin(a));

    % 2. Solve for Trim (Cm = 0)
    alpha_interval = [deg2rad(0.1), deg2rad(45)];

    try
        alpha_trim = fzero(Cm_fun, alpha_interval);
    catch
        % If fzero fails, no trim point exists in that range
        alpha_trim = NaN;
    end

       % 4. Calculate Flight Speed and Thrust at Trim
    V_trim = NaN;
    Thrust_req = NaN;

    if ~isnan(alpha_trim)
        CL_t = C_L_interp(alpha_trim);
        CD_t = C_D_interp(alpha_trim);
        
        if CL_t > 0
            % Normal (CN) and Axial (CA) coefficients
            CN_t = CL_t * cos(alpha_trim) + CD_t * sin(alpha_trim);
            CA_t = -CL_t * sin(alpha_trim) + CD_t * cos(alpha_trim);
            
            % Speed to balance weight component normal to the body
            V_trim = sqrt((2 * W * cos(alpha_trim)) / (rho * S_ref * CN_t));
            
            % Thrust to balance aerodynamic axial drag + gravity component
            q_t = 0.5 * rho * V_trim^2;
            Thrust_req = (q_t * S_ref * CA_t) + (W * sin(alpha_trim));
        end
    end
end
