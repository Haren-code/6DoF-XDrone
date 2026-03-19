%% --- INTEGRATION WITH LIVE PLOT ---
function options_integration = options_premaker(live_plotting_calculation, sim)

    if live_plotting_calculation == false
        % Simple plotting without residuals (faster)
        options_integration = odeset('RelTol',1e-6,'AbsTol',1e-9,'Events', @hit_ground);
    else
        % Plotting with residuals of the x states (slow)
        fig = figure('Name','Live dx evolution','NumberTitle','off');
        tiledlayout(fig, 4, 4, 'TileSpacing','tight'); 
        sgtitle('Time evolution of dx components')
        
        % Store parameters & plot handles for access inside OutputFcn
        setappdata(fig, 'sim', sim);
        %setappdata(fig, 'x0', x0);
        
        % Define the names of your 13 derivatives
        dx_names = {'du (Axial Accel)', 'dv (Lat Accel)', 'dw (Norm Accel)', ...
                    'dp (Roll Accel)', 'dq_ang (Pitch Accel)', 'dr (Yaw Accel)', ...
                    'dq0_rate', 'dq1_rate', 'dq2_rate', 'dq3_rate', ...
                    'dx (Inertial Vel X)', 'dy (Inertial Vel Y)', 'dz (Inertial Vel Z)'};

        % Create line objects using animatedline for speed
        for i = 1:13
            ax(i) = nexttile;
            title(ax(i), dx_names{i});
            xlabel(ax(i), 't [s]');
            ylabel(ax(i), 'Rate of Change');
            grid(ax(i), 'on');
            
            % animatedline is specifically designed for streaming live data
            h(i) = animatedline(ax(i), 'Color', 'b', 'Marker', '.', 'LineStyle', 'none');
        end
        
        setappdata(fig, 'h', h);
        setappdata(fig, 'ax', ax);
        
        options_integration = odeset('RelTol',1e-6,'AbsTol',1e-9,'Events',@hit_ground,...
            'OutputFcn', @(t,y,flag) live_dx_plot(t,y,flag,fig));
    end
end