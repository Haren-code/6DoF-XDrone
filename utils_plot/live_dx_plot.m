function status = live_dx_plot(t, y, flag, fig)
% Live plotting function for dx components during ODE45 integration
    
    status = 0; % always continue integration
    
    if strcmp(flag, 'init')
        % nothing special at init
    elseif isempty(flag)
        sim = getappdata(fig, 'sim');
        h = getappdata(fig, 'h');
    
        for i = 1:length(t)
            % Re-evaluate physics to get the derivatives (dx)
            dx = sixDoF_wrapper(t(i), y(:,i), sim); 
            
            % Instantly add points to the animated lines
            for j = 1:min(length(dx), numel(h))
                addpoints(h(j), t(i), dx(j));
            end
        end
        drawnow limitrate nocallbacks
    elseif strcmp(flag, 'done')
        disp('Integration finished.');
    end
end