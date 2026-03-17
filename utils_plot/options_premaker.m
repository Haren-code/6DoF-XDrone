    %% --- INTEGRATION WITH LIVE PLOT ---
function options_integration = options_premaker(live_plotting_residuals, sim)

    if live_plotting_residuals == false
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
        
        % Create line objects for speed (no replotting overhead)
        for i = 1:13
            ax(i) = nexttile;
            hold(ax(i), 'on');
            title(ax(i), sprintf('dx(%d)', i));
            xlabel(ax(i), 't [s]');
            ylabel(ax(i), sprintf('dx_%d', i));
            grid(ax(i), 'on');
            h(i) = plot(ax(i), NaN, NaN, 'b.');
        end
        
        setappdata(fig, 'h', h);
        setappdata(fig, 'ax', ax);
        
        options_integration = odeset('RelTol',1e-6,'AbsTol',1e-9,'Events',@hit_ground,...
            'OutputFcn', @(t,y,flag) live_dx_plot(t,y,flag,fig));
    end
end