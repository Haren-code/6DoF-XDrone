function xZylo_visualization
    % --- Initial Setup ---
    clear
    clc
    close all
    xzylo_wanted = 0;

    if xzylo_wanted == 1
        [aero,~] = xzylo();
      
        % --- 1. PARAMETERS SETUP ---
        
        % Define the cylinder's physical properties
        R_ext = aero.b;    % Exterior Radius (m)
        L = aero.c;        % Length of the cylinder (m)
        Weight = aero.m*aero.g;   % Weight for gravity calculation (N)
        force_scale = 10;
        
        % COG position (fixed at COG_x = 1 in the inertial frame)
        COG_x = aero.CoG(1);

    else 
    
        % Define the cylinder's physical properties
        R_ext = 1;    % Exterior Radius (m)
        L = 1;        % Length of the cylinder (m)
        Weight = 3;   % Weight for gravity calculation (N)
        force_scale = 10;
        
        % COG position (fixed at COG_x = 1 in the inertial frame)
        COG_x = 0.3;
    end


    % --- 2. FIGURE AND AXES SETUP (Maximized, No Box, Vis3d) ---
    
    % Create the main figure window, setting it to maximized
    h_fig = figure('Name', 'Cylinder Rotation and Force Visualization', ...
                   'NumberTitle', 'off', ...
                   'WindowState', 'maximized'); 
    
    % Create the 3D axes
    h_ax = axes('Parent', h_fig); 
    hold(h_ax, 'on');
    axis(h_ax, 'equal');
    
    % --- AXES APPEARANCE COMMANDS ---
    box(h_ax, 'off');  % Removes the outer perimeter frame
    grid(h_ax, 'on');  % Keeps the grid visible
    axis(h_ax, 'vis3d'); % Preserves aspect ratio and prevents clipping
    
    xlabel(h_ax, 'X (Flight Direction)');
    ylabel(h_ax, 'Y (Side)');
    zlabel(h_ax, 'Z (Up)');
    title(h_ax, 'Cylinder with Aerodynamic Forces');
    
    view(h_ax, 45, 30); % Set an initial 3D view
    
    % Set axis limits
    ax_lim = L * 3;
    xlim(h_ax, [-ax_lim, ax_lim]);
    ylim(h_ax, [-ax_lim, ax_lim]);
    zlim(h_ax, [-ax_lim, ax_lim]);
    
    % --- 3. CYLINDER GEOMETRY DEFINITION ---
    
    % Generate coordinates for a cylinder along the x-axis (from -L to 0)
    [Z_ext, Y_ext, X_ext] = cylinder(R_ext, 50);
    X_ext = COG_x - X_ext * L; 
    
    % Create the cylinder surface object 
    h_cyl = surf(h_ax, X_ext, Y_ext, Z_ext, 'FaceColor', [0.7 0.7 0.7], ...
                 'EdgeColor', 'k', 'FaceAlpha', 0.6);

    % Define the cylinder's central axis in its local frame (X-axis from -L/2 to L/2)
    Axis_X = [-1.5*L, 0.5*L] + COG_x ;
    Axis_Y = [0, 0];
    Axis_Z = [0, 0];
    
    % 1. Dotted Red Axis
    h_axis = plot3(h_ax, Axis_X, Axis_Y, Axis_Z, ...
                   'r--', 'LineWidth', 1.5);
               
    % 2. Thick Red COG Point
    % The COG is defined at X = COG_x (1 meter from the origin in the local frame)
    h_cog = plot3(h_ax, 0, 0, 0, ...
                  'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    
    % --- 4. SLIDERS SETUP (UI CONTROLS) - Bottom-Left Corner ---
    
    min_angle = -90;
    max_angle = 90;
    initial_alpha = 0;
    initial_beta = 0;
    
    % Define layout parameters (in pixels)
    slider_width  = 300;
    slider_height = 20;
    label_width   = 150;
    label_height  = 20;
    value_width   = 60; 
    
    % Bottom-left corner layout position
    base_y = 50;  
    base_x = 30;  
    
    % --- Alpha Slider (Pitch - rotation around Y-axis) ---
    uicontrol('Style', 'text', ...
              'String', 'Pitch (α) [deg]:', ...
              'Position', [base_x, base_y + 30, label_width, label_height], ...
              'HorizontalAlignment', 'left');
    
    h_slider_alpha = uicontrol('Style', 'slider', ...
                               'Min', min_angle, 'Max', max_angle, 'Value', initial_alpha, ...
                               'Position', [label_width - 30, base_y + 30, slider_width, slider_height], ...
                               'Callback', @updateVisualization);
    
    % Direct Input Box for Alpha 
    h_edit_alpha = uicontrol('Style', 'edit', ...
                             'String', sprintf('%.1f', initial_alpha), ...
                             'Position', [label_width + slider_width - 20, base_y + 30, value_width, label_height], ...
                             'Callback', @updateEditAlpha); 
    
    % --- Beta Slider (Yaw - rotation around Z-axis) ---
    uicontrol('Style', 'text', ...
              'String', 'Yaw (β) [deg]:', ...
              'Position', [base_x, base_y, label_width, label_height], ...
              'HorizontalAlignment', 'left');
    
    h_slider_beta = uicontrol('Style', 'slider', ...
                              'Min', min_angle, 'Max', max_angle, 'Value', initial_beta, ...
                              'Position', [label_width - 30, base_y, slider_width, slider_height], ...
                              'Callback', @updateVisualization);
    
    % Direct Input Box for Beta
    h_edit_beta = uicontrol('Style', 'edit', ...
                            'String', sprintf('%.1f', initial_beta), ...
                            'Position', [label_width + slider_width - 20, base_y, value_width, label_height], ...
                            'Callback', @updateEditBeta); 
    
    % --- Reset Button ---
    uicontrol('Style', 'pushbutton', ...
              'String', 'Reset (0°)', ...
              'Position', [base_x, base_y - 40, 100, 30], ...
              'Callback', @resetAngles); 
             
    % --- 5. INITIAL FORCE ARROWS ---
    
    % Gravity (Green) 
    h_grav = quiver3(h_ax, 0, 0, 0, 0, 0, -Weight/force_scale, 'g', 'LineWidth', 2, 'MaxHeadSize', 0.5,'AutoScale', 'off');
    % Lift (Red)
    h_lift = quiver3(h_ax, 0, 0, 0, 0, 0, 0, 'r', 'LineWidth', 2, 'MaxHeadSize', 0.5,'AutoScale', 'off');
    % Side Force (Orange) 
    h_side = quiver3(h_ax, 0, 0, 0, 0, 0, 0, 0, 'Color', [1 0.6 0], 'LineWidth', 2, 'MaxHeadSize', 0.5, 'AutoScale', 'off');
    % Total Lift/Side (Magenta)
    h_total_lat = quiver3(h_ax, 0, 0, 0, 0, 0, 0, 0, 'm', 'LineWidth', 2, 'MaxHeadSize', 0.5,'AutoScale', 'off');
    % Drag (Cyan/Light Blue) 
    h_drag = quiver3(h_ax, 0, 0, 0, 0, 0, 0, 0, 'c', 'LineWidth', 2, 'MaxHeadSize', 0.5,'AutoScale', 'off');
    
       
    % --- 6. FIXED WIND ARROWS (Wall of Blue Arrows) ---
    
    R_grid = 2*L; 
    N_arrows = 5; 
    
    Y_positions = linspace(-R_grid, R_grid, N_arrows);
    Z_positions = linspace(-R_grid, R_grid, N_arrows);
    
    wind_start_x = 3 * L; 
    wind_end_x = 2.5 * L;   
    wind_vector_x = wind_end_x - wind_start_x; 
    
    for i = 1:N_arrows
        for j = 1:N_arrows
            y_start = Y_positions(i);
            z_start = Z_positions(j);
            
            % Draw the arrow using quiver3
            quiver3(h_ax, wind_start_x, y_start, z_start, ...
                    wind_vector_x, 0, 0, 0, ... 
                    'b', 'LineWidth', 1.5, 'MaxHeadSize', 0.5);
        end
    end
    
    % Add legend for forces
    legend([h_grav, h_lift, h_side, h_total_lat, h_drag, h_cyl], ...
           {'Gravity (Const.)', 'Lift (f(\alpha))', 'Side (f(\beta))', ...
            'Total Lateral', 'Drag', 'Cylinder'}, ...
           'Location', 'best');
    
    % Call the update function once to draw the initial state
    updateVisualization();
  

    
    % --- 7. UPDATE FUNCTION (Main Visualization Logic) ---
    
    function updateVisualization(~, ~)
        % Get current angles from sliders (in degrees)
        alpha_deg = get(h_slider_alpha, 'Value'); 
        beta_deg = get(h_slider_beta, 'Value');   
        
        % Update the edit box strings to synchronize with the slider
        % This ensures the text box shows the current slider value
        set(h_edit_alpha, 'String', sprintf('%.1f', alpha_deg));
        set(h_edit_beta, 'String', sprintf('%.1f', beta_deg));
        
        % Convert angles to radians
        alpha_rad = deg2rad(-alpha_deg);
        beta_rad = deg2rad(beta_deg);
        
        % --- A. CYLINDER ROTATION ---
        
        % Rotation Matrices:
        Ry = [cos(alpha_rad), 0, sin(alpha_rad);
              0, 1, 0;
              -sin(alpha_rad), 0, cos(alpha_rad)];
        
        Rz = [cos(beta_rad), -sin(beta_rad), 0;
              sin(beta_rad), cos(beta_rad), 0;
              0, 0, 1];
        
        R_total = Rz * Ry;
        
        % Flatten the cylinder coordinates
        coords_cyl = [X_ext(:)'; Y_ext(:)'; Z_ext(:)'];
        
        % Apply rotation
        coords_rotated = R_total * coords_cyl;
        
        % Reshape the rotated coordinates back
        X_rot = reshape(coords_rotated(1, :), size(X_ext));
        Y_rot = reshape(coords_rotated(2, :), size(Y_ext));
        Z_rot = reshape(coords_rotated(3, :), size(Z_ext));
        
        % Update the cylinder's surface data
        set(h_cyl, 'XData', X_rot, 'YData', Y_rot, 'ZData', Z_rot);
      
        % 2. Body Axis and COG Point (NEW TRANSFORMATION)
        
        % Body Axis end points in Local Frame (X: -L/2 to L/2, Y/Z: 0)
        Axis_coords_body = [Axis_X; Axis_Y; Axis_Z];
        
        % Apply rotation to the axis
        Axis_coords_inertial = R_total * Axis_coords_body;
        
        % Update the Axis
        set(h_axis, 'XData', Axis_coords_inertial(1, :) , ...
                    'YData', Axis_coords_inertial(2, :) , ...
                    'ZData', Axis_coords_inertial(3, :) );
                
        % COG Position in Local Frame (X: COG_x, Y/Z: 0)
        COG_coords_body = [0; 0; 0];
        
        % Apply rotation to the COG
        COG_coords_inertial = R_total * COG_coords_body;
        
        % Update the COG Marker
        set(h_cog, 'XData', COG_coords_inertial(1, 1) , ...
                   'YData', COG_coords_inertial(2, 1) , ...
                   'ZData', COG_coords_inertial(3, 1) );

        % Calculate distance from (0,0,0) to the new COG coordinates
        X_abs_cog = COG_coords_inertial(1, 1) ;
        Y_abs_cog = COG_coords_inertial(2, 1) ;
        Z_abs_cog = COG_coords_inertial(3, 1) ;
        Distance_to_Origin = norm([X_abs_cog, Y_abs_cog, Z_abs_cog]); % norm is sqrt(x^2 + y^2 + z^2)
        
        % Print the required message to the Command Window
        fprintf('The distance between the COG and origin is %.4f m, where the COG coordinates are [%.4f, %.4f, %.4f]\n', ...
                Distance_to_Origin, X_abs_cog, Y_abs_cog, Z_abs_cog);
            
        
        % --- B. FORCE CALCULATION (In the BODY Frame) ---
        
        % Gravity is constant in the inertial frame
        G_inertial = [0, 0, -Weight]; 
        set(h_grav, 'UData', G_inertial(1)/force_scale, 'VData', G_inertial(2)/force_scale, 'WData', G_inertial(3)/force_scale);

        % Aerodynamic Forces (Defined in the BODY Frame)

        num = sqrt( sin(beta_rad).^2 + (sin(alpha_rad).^2).*(cos(beta_rad).^2) );
        den = cos(alpha_rad).*cos(beta_rad);
        theta_total_rad = atan2(num, den);
        theta_total_deg = rad2deg(theta_total_rad);

        fprintf('The total angle of attack is: %.4f \n', theta_total_deg);

        if xzylo_wanted == 1
            Lift_mag = aero.CL(alpha_rad,beta_rad);
            Side_force_mag = aero.CY(alpha_rad,beta_rad);          
            Drag_mag = aero.CD(theta_total_rad,beta_rad);
        else
            Lift_mag = -0.25*alpha_deg;
            Side_force_mag = 0.25*beta_deg;
            Drag_mag = 0.4*theta_total_deg;
        end

        F_Lift = [0, 0, -Lift_mag];
        F_Side = [0, Side_force_mag, 0];
        F_Drag = [-Drag_mag, 0, 0];
        F_TotalLat = F_Lift + F_Side;

        % --- D. UPDATE FORCE ARROWS ---

        set(h_lift, 'UData', F_Lift(1)/force_scale, 'VData', F_Lift(2)/force_scale, ...
                    'WData', F_Lift(3)/force_scale);

        set(h_side, 'UData', F_Side(1)/force_scale, 'VData', F_Side(2)/force_scale, ...
                    'WData', F_Side(3)/force_scale);

        set(h_total_lat, 'UData', F_TotalLat(1)/force_scale, 'VData', F_TotalLat(2)/force_scale, ...
                         'WData', F_TotalLat(3)/force_scale);

        set(h_drag, 'UData', F_Drag(1)/force_scale, 'VData', F_Drag(2)/force_scale, ...
                    'WData', F_Drag(3)/force_scale);
    end

    % --- FUNCTIONS for direct text input ---
    
    function updateEditAlpha(~, ~)
        % Read value from edit box
        new_alpha = str2double(get(h_edit_alpha, 'String'));
        
        % Validate value against limits
        if isnan(new_alpha) || new_alpha < min_angle || new_alpha > max_angle
            % If input is invalid, reset edit box and slider to current slider value
            set(h_edit_alpha, 'String', sprintf('%.1f', get(h_slider_alpha, 'Value')));
            return;
        end
        
        % Update slider value and call the main visualization function
        set(h_slider_alpha, 'Value', new_alpha);
        updateVisualization();
    end

    function updateEditBeta(~, ~)
        % Read value from edit box
        new_beta = str2double(get(h_edit_beta, 'String'));
        
        % Validate value against limits
        if isnan(new_beta) || new_beta < min_angle || new_beta > max_angle
            % If input is invalid, reset edit box and slider to current slider value
            set(h_edit_beta, 'String', sprintf('%.1f', get(h_slider_beta, 'Value')));
            return;
        end
        
        % Update slider value and call the main visualization function
        set(h_slider_beta, 'Value', new_beta);
        updateVisualization();
    end
    
    % --- RESET FUNCTION ---
    function resetAngles(~, ~)

        set(h_slider_alpha, 'Value', 0.0);
        set(h_slider_beta, 'Value', 0.0);
        
        updateVisualization();

    end
end