function plot_results(t, x, sim, SIM_DATA)
    
    % --- Handle default params ---
    if nargin < 3; sim.options = struct(); end
    if ~isfield(sim.options, 'body_plotting'); sim.options.body_plotting = false; end
    if ~isfield(sim.options, 'test_plotting'); sim.options.test_plotting = false; end
    if ~isfield(sim.options, 'non_rotating_plotting'); sim.options.non_rotating_plotting = false; end

    % ---------------------- Extract states ------------------------------
    u = x(:,1); v = x(:,2); w = x(:,3);         % body velocities
    p = x(:,4); q = x(:,5); r = x(:,6);         % angular rates
    q_quat = x(:,7:10);                         % quaternions
    pos = x(:,11:13);                           % inertial positions
        
    % ------ Convert quaternions to Euler angles (ZYX) --------
    eul = quat2eul(q_quat, 'ZYX');              
    psi = eul(:,1); theta = eul(:,2); phi = eul(:,3);
    
    %% ==== CALCULATION: NON-ROTATING BODY FRAME (De-spun) ====
    % We apply the inverse roll rotation to vectors in the body frame
    
    % Velocities (u is x-axis, so u_nr = u)
    u_nr = u;
    v_nr = v .* cos(phi) - w .* sin(phi);
    w_nr = v .* sin(phi) + w .* cos(phi);
    alpha_nr = atan2(w_nr, u_nr);   % AoA
    V_total = sqrt(u_nr.^2 + w_nr.^2 + v_nr.^2);
    beta_nr = asin(v_nr ./ V_total);
    q_nr = q .* cos(phi) - r .* sin(phi);
    r_nr = q .* sin(phi) + r .* cos(phi);

    %% ====================== STANDARD PLOTS =============================
    % These plots show the fundamental state of the vehicle
    
    % 1. Body linear velocities (States)
    figure('Name','Body Velocities ');
    plot(t,u,'r',t,v,'g',t,w,'b','LineWidth',1.5);
    xlabel('Time [s]'); ylabel('Body velocities [m/s]');
    legend('u','v','w'); grid on; title('Body-frame linear velocities');

    % 2. Angular velocities (States)
    figure('Name','Angular Velocities');
    plot(t,p/(2*pi),'r',t,q/(2*pi),'g',t,r/(2*pi),'b','LineWidth',1.5);
    xlabel('Time [s]'); ylabel('Angular rates [Hz]');
    legend('p','q','r'); grid on; title('Body angular rates');

    % 3. Euler angles
    figure('Name','Euler Angles (without roll)');
    plot(t, rad2deg([theta psi]), 'LineWidth', 1.5);
    xlabel('Time [s]'); ylabel('Euler angles [deg]');
    legend('\theta (pitch)','\psi (yaw)');
    grid on; title('Euler angles');

    % 4. Aerodynamic Angles
    figure('Name','Aerodynamic Angles')
    hold on
    % plot(SIM_DATA.t, rad2deg(SIM_DATA.alpha), 'LineWidth', 1.5);
    % plot(SIM_DATA.t, rad2deg(SIM_DATA.beta), 'LineWidth', 1.5);
    plot(SIM_DATA.t, rad2deg(SIM_DATA.alpha_equivalent), 'LineWidth', 1.5);
    xlabel('Time [s]'); ylabel('Aerodynamic angles [deg]');
    legend('\alpha_{equivalent}');
    grid on; title('Aerodynamic angles');
    
    % 5. Aero coefficients
    figure('Name', 'Aero coefficients')
    hold on
    plot(SIM_DATA.t, SIM_DATA.CL, 'LineWidth', 1.5);
    plot(SIM_DATA.t, SIM_DATA.CD, 'LineWidth', 1.5);
    plot(SIM_DATA.t, SIM_DATA.Cm, 'LineWidth', 1.5);
    xlabel('Time [s]'); ylabel('Aerodynamic coefficients');
    legend('C_L','C_D','C_m');
    grid on; title('Aerodynamic coefficients (around CoG)');

    % 6. CoP vs CoG
    figure('Name', 'CoP position')
    hold on
    plot(SIM_DATA.t, SIM_DATA.CoP_fraction * 100, 'LineWidth', 1.5);
    yline(sim.prop.Xzylo.percentage_CoG(1) * 100, '-', 'LineWidth', 1.5);
    xlabel('Time [s]'); ylabel('Position from LE in percentage');
    legend('CoP position','CoG position');
    grid on; title('CoP position variation');

    % 7. General Aero forces (Wind/Stability axes)
    figure('Name','Aero forces')
    hold on
    plot(SIM_DATA.t, SIM_DATA.Lift_sdslip, 'LineWidth', 1.5)
    plot(SIM_DATA.t, SIM_DATA.Drag_sdslip, 'LineWidth', 1.5)
    plot(SIM_DATA.t, SIM_DATA.Normal_force, 'LineWidth', 1.5)
    plot(SIM_DATA.t, SIM_DATA.Axial_force, 'LineWidth', 1.5)
    xlabel('Time [s]'); ylabel('Aero forces');
    legend('Lift','Drag','Normal','Axial');
    grid on; title('Aero forces in wind frame and no-sideslip frames');

    % 8. 3D Trajectory
    figure('Name','Trajectory');
    plot3(pos(:,1), pos(:,2), pos(:,3), 'r','LineWidth',1.5);
    xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
    grid on; axis equal; title('Trajectory in inertial frame');
    set(gca, 'YDir', 'reverse', 'ZDir', 'reverse');  % flip East and Down visually
    

    %% ====================== NON-ROTATING FRAME PLOTS ===================
    if sim.options.non_rotating_plotting
        
        % 1. Non-Rotating Velocities
        figure('Name','Non-Rotating Velocities');
        plot(t, u_nr, 'r', t, v_nr, 'g', t, w_nr, 'b', 'LineWidth', 1.5);
        xlabel('Time [s]'); ylabel('Velocities [m/s]');
        legend('u_{nr}', 'v_{nr}', 'w_{nr}'); 
        grid on; title('Velocities in Non-Rotating Frame (De-spun)');
        
        % 2. Non-Rotating Forces
        figure('Name','Non-Rotating Forces');
        hold on;
        plot(SIM_DATA.t, SIM_DATA.Fx_nr, 'LineWidth', 1.5);
        plot(SIM_DATA.t, SIM_DATA.Fy_nr, 'LineWidth', 1.5);
        plot(SIM_DATA.t, SIM_DATA.Fz_nr, 'LineWidth', 1.5);
        xlabel('Time [s]'); ylabel('Forces [N]');
        legend('F_{x,nr}', 'F_{y,nr}', 'F_{z,nr}');
        grid on; title('Forces in Non-Rotating Frame (De-spun)');
        
        % 3. Non-Rotating Moments
        figure('Name','Non-Rotating Moments');
        hold on;
        plot(SIM_DATA.t, SIM_DATA.Mx_nr, 'LineWidth', 1.5);
        plot(SIM_DATA.t, SIM_DATA.My_nr, 'LineWidth', 1.5);
        plot(SIM_DATA.t, SIM_DATA.Mz_nr, 'LineWidth', 1.5);
        xlabel('Time [s]'); ylabel('Moments [Nm]');
        legend('M_{x,nr}', 'M_{y,nr}', 'M_{z,nr}');
        grid on; title('Moments in Non-Rotating Frame (De-spun)');

        % 4. Aerodynamic Angles
        figure('Name','Non-Rotating Angles');
        plot(t, rad2deg(alpha_nr), 'r', t, rad2deg(beta_nr), 'g','LineWidth', 1.5);
        xlabel('Time [s]'); ylabel('Angles');
        legend('alpha_{nr}', 'beta_{nr}'); 
        grid on; title('Aero Angles in Non-Rotating Frame (De-spun)');

        % 5.Non-Rotating Angular Velocities
        figure('Name','Non-Rotating Angular Velocities');
        plot(t, p/(2*pi), 'r', t, q_nr/(2*pi), 'g', t, r_nr/(2*pi), 'b', 'LineWidth', 1.5);
        xlabel('Time [s]'); ylabel('Angular Velocities [Hz]');
        legend('p_{nr}', 'q_{nr}', 'r_{nr}'); 
        grid on; title('Angular Velocities in Non-Rotating Frame (De-spun)');
    end

    %% ====================== BODY FRAME PLOTS ===========================
    % Triggered if body_plotting = true
    if sim.options.body_plotting
        
        % Body Aerodynamic Forces
        figure('Name','Body aerodynamic forces')
        hold on
        plot(SIM_DATA.t, SIM_DATA.Fx_aero_b, 'LineWidth', 1.5)
        plot(SIM_DATA.t, SIM_DATA.Fy_aero_b, 'LineWidth', 1.5)
        plot(SIM_DATA.t, SIM_DATA.Fz_aero_b, 'LineWidth', 1.5)
        xlabel('Time [s]'); ylabel('Body aerodynamic forces');
        legend('F_{aero,x}','F_{aero,y}','F_{aero,z}');
        grid on; title('Aero forces in body frame');
        
        % Total Body Forces (inc Gravity)
        figure('Name','Total body forces')
        hold on
        plot(SIM_DATA.t, SIM_DATA.Fx_b, 'LineWidth', 1.5)
        plot(SIM_DATA.t, SIM_DATA.Fy_b, 'LineWidth', 1.5)
        plot(SIM_DATA.t, SIM_DATA.Fz_b, 'LineWidth', 1.5)
        xlabel('Time [s]'); ylabel('Total body forces');
        legend('F_{body,x}','F_{body,y}','F_{body,z}');
        grid on; title('Body forces (including gravity)');
        
        % Body Moments
        figure('Name', 'Body moments')
        hold on
        plot(SIM_DATA.t, SIM_DATA.Mx_b, 'LineWidth', 1.5);
        plot(SIM_DATA.t, SIM_DATA.My_b, 'LineWidth', 1.5);
        plot(SIM_DATA.t, SIM_DATA.Mz_b, 'LineWidth', 1.5);
        xlabel('Time [s]'); ylabel('Body moments');
        legend('M_x','M_y','M_z');
        grid on; title('Moments in body frame');
    end


    %% ====================== TEST / VERIFICATION PLOTS ==================
    % Triggered if test_plotting = true
    if sim.options.test_plotting
        
        % Test: No Sideslip Velocities
        figure('Name', 'No sideslip velocities')
        hold on
        plot(SIM_DATA.t, SIM_DATA.u_no_sdslip, 'LineWidth', 1.5);
        plot(SIM_DATA.t, SIM_DATA.v_no_sdslip, 'LineWidth', 1.5);
        plot(SIM_DATA.t, SIM_DATA.w_no_sdslip, 'LineWidth', 1.5);
        xlabel('Time [s]'); ylabel('Body velocities in no-sideslip frame [m/s]');
        legend('u_{no,sdslip}','v_{no,sdslip}','w_{no,sdslip}');
        grid on; title('Test if v = 0');
        
        % Test: Velocity Norm Calculation
        figure('Name', 'Test velocity norm')
        hold on
        plot(SIM_DATA.t, SIM_DATA.V_no_sdslip, 'LineWidth', 1.5);
        plot(SIM_DATA.t, SIM_DATA.V, 'LineWidth', 1.5);
        plot(t, sqrt(u.^2 + v.^2 + w.^2), 'LineWidth', 1.5);
        xlabel('Time [s]'); ylabel('Velocity Norm');
        legend('V_{sideslip}','V_{norm}','V_{sqrt,calculated}');
        grid on; title('Velocity norm verification');
        
        % Test: No-sideslip frame moments (Verification)
        figure('Name', 'No-sideslip frame moments')
        hold on
        plot(SIM_DATA.t, SIM_DATA.Mx_no_sdslip, 'LineWidth', 1.5);
        plot(SIM_DATA.t, SIM_DATA.My_no_sdslip, 'LineWidth', 1.5);
        plot(SIM_DATA.t, SIM_DATA.Mz_no_sdslip, 'LineWidth', 1.5);
        xlabel('Time [s]'); ylabel('Moments in no-sideslip frame');
        legend('M_{sdslip,x}','M_{sdslip,y}','M_{sdslip,z}');
        grid on; title('No-sideslip moments (Mx and Mz should be zero)');
    end
end