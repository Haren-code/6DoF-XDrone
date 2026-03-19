clear; clc; close all;

%% ========================================================================
%%                    PART 1: GEOMETRY & CHORD SCALING
%% ========================================================================

%% ============ USER INPUTS ============
writetable_value = true; 

% -- Design Conditions --
V_inf_design = 15;      % [m/s]
RPM_design   = 1500;    % [RPM]
Omega_design = RPM_design*2*pi/60;

% -- Geometry Dimensions --
R_hub  = 11e-3;         % [m]
R_duct = 56.7e-3;       % [m]
hub_h  = 0.03;          % [m]
duct_h = 0.115;         % [m]
c_base = hub_h - 2e-3;  % [m] Base Chord (before scaling)
B  = 2;                 % Number of Blades

% --- BLADE MODE ---
% 1 = Constant Z (Chord Projection) -> Chord INCREASES with twist
% 2 = Constant Chord (Twisted)      -> Chord is constant
blade_mode = 1; 

rotation_position = 0;  % 0 = LE (Since LE is at 0,0)
alpha_d = deg2rad(0);   % ZERO LIFT DESIGN
N_span  = 50;           % Grid resolution
rho     = 1.225;        % Air Density

%% =====================================
x_span = linspace(R_hub, R_duct, N_span);

% --- TWIST CALCULATION ---
% Phi = Angle from Vertical (Axial) axis.
phi_design = atan((Omega_design * x_span) ./ V_inf_design);
% Beta = Physical Twist from Vertical
beta_phys = phi_design + alpha_d; 

%% -------- IMPORT GEOMETRY (NACA0012) --------
airfoil_filename = 'NACA0012.dat';

if isfile(airfoil_filename)
    raw_data = readmatrix(airfoil_filename, 'NumHeaderLines', 2);
    split_idx = find(diff(raw_data(:,1)) < -0.5, 1);
    
    if isempty(split_idx)
        error('Could not identify Upper/Lower split in airfoil file.');
    else
        Upper_X = raw_data(1:split_idx, 1); Upper_Y = raw_data(1:split_idx, 2);
        Lower_X = raw_data(split_idx+1:end, 1); Lower_Y = raw_data(split_idx+1:end, 2);
        X_segs = [Upper_X; 1.0; flipud(Lower_X)];
        Y_segs = [Upper_Y; 0.0; flipud(Lower_Y)];
        X_norm = X_segs'; Y_norm = Y_segs';
        idx_LE = 1; idx_TE = length(Upper_X) + 1; 
    end
    fprintf('Loaded Airfoil. LE Index: %d, TE Index: %d\n', idx_LE, idx_TE);
else
    error('Airfoil file %s not found!', airfoil_filename);
end

% Base Airfoil (Unscaled, Unrotated)
X_air_base = (X_norm - rotation_position) * c_base; 
Y_air_base = Y_norm * c_base;

%% -------- CALCULATE SCALED CHORD DISTRIBUTION --------
fprintf('Calculating Chord Scaling based on Mode %d...\n', blade_mode);
c_dist = zeros(size(x_span)); 
beta_root = beta_phys(1) + pi;
Z_root = X_air_base .* cos(beta_root) + Y_air_base .* sin(beta_root);
H_ref = abs(Z_root(idx_TE) - Z_root(idx_LE));

for i = 1:N_span
    beta = beta_phys(i) + pi;
    Z_curr = X_air_base .* cos(beta) + Y_air_base .* sin(beta);
    H_curr = abs(Z_curr(idx_TE) - Z_curr(idx_LE));
    if blade_mode == 1 
        if H_curr > 1e-6; scale_factor = H_ref / H_curr; else; scale_factor = 1; end
    else
        scale_factor = 1;
    end
    c_dist(i) = c_base * scale_factor;
end

%% ========================================================================
%%            PART 2: INDUCTION ANALYSIS (TURBINE BEM)
%% ========================================================================
fprintf('\n--- Running Induction Check (Turbine BEM) ---\n');

% Load Polar Data
csv_filename = 'xf-n0012-il-200000.csv';
PolarData = readtable(csv_filename, 'NumHeaderLines', 10);
[Alpha_data, uI] = unique(PolarData.Alpha);
Cl_data = PolarData.Cl(uI); Cd_data = PolarData.Cd(uI);

a_dist = zeros(1,N_span-1); ap_dist = zeros(1,N_span-1); 
alpha_check = zeros(1,N_span-1); r_plot = zeros(1,N_span-1);

for k = 1:N_span-1
    r_local = (x_span(k) + x_span(k+1))/2; r_plot(k) = r_local;
    c_local = (c_dist(k) + c_dist(k+1))/2;
    beta_local = (beta_phys(k) + beta_phys(k+1))/2;
    sigma = (B * c_local) / (2 * pi * r_local);
    
    a = 0; ap = 0; alpha_deg = 0; 
    
    if k == 1
        a_dist(k) = NaN; ap_dist(k) = NaN; alpha_check(k) = NaN;
        continue; 
    end
    
    if r_local > (R_hub + 1e-4)
        for iter = 1:100 
            denom = V_inf_design * (1 - a); 
            if abs(denom)<1e-4; denom=1e-4; end 
            
            phi = atan( (Omega_design * r_local * (1 + ap)) / denom );
            sin_p = sin(phi); cos_p = cos(phi);
            if abs(sin_p) < 1e-3; sin_p = 1e-3; end
            
            % Prandtl Loss
            f_tip = (B/2) * (R_duct - r_local) / (r_local * abs(sin_p));
            F_tip = (2/pi) * acos(min(1, max(-1, exp(-f_tip))));
            f_hub = (B/2) * (r_local - R_hub) / (r_local * abs(sin_p));
            F_hub = (2/pi) * acos(min(1, max(-1, exp(-f_hub))));
            F = F_tip * F_hub; 
            if F < 1e-4; F = 1e-4; end 
            
            alpha_deg = rad2deg(beta_local - phi);
            
            Cl = interp1(Alpha_data, Cl_data, alpha_deg, 'linear', 'extrap');
            Cd = interp1(Alpha_data, Cd_data, alpha_deg, 'linear', 'extrap');
            
            % TURBINE PHYSICS: Swapped coefficients
            Cn = Cl * sin_p + Cd * cos_p; % Axial 
            Ct = Cl * cos_p - Cd * sin_p; % Tangential
            
            if abs(Cn) < 1e-5
                a_new = 0; ap_new = 0;
            else
                CT_local = sigma * Cn * (1-a)^2 / (sin_p^2);
                % Standard BEM induction calculations
                denom_a = (4 * F * sin_p^2) / (sigma * Cn);
                a_new = 1 / (denom_a + 1);
                
                denom_ap = (4 * F * sin_p * cos_p) / (sigma * Ct);
                ap_new = 1 / (denom_ap - 1);
            end
            
            relax = 0.1;
            a = (1-relax)*a + relax*a_new; 
            ap = (1-relax)*ap + relax*ap_new;
        end
    end
    
    if a > 1; a = 1; end; if a < -0.5; a = -0.5; end
    a_dist(k) = a; ap_dist(k) = ap; alpha_check(k) = alpha_deg;
end

a_dist(1) = a_dist(2);
ap_dist(1) = ap_dist(2);
alpha_check(1) = alpha_check(2);

%% ============== FIGURE: INDUCTION ANALYSIS ==================
figure('Color','w','Name','Design Point Check','WindowStyle','docked');
subplot(2,1,1); plot(r_plot, a_dist, 'r-o', 'LineWidth', 1.5); grid on; title('Axial Induction a'); xlabel('Span [m]');
subplot(2,1,2); plot(r_plot, alpha_check, 'k-d', 'LineWidth', 1.5); yline(rad2deg(alpha_d),'r--','Target'); grid on; xlabel('Span [m]'); ylabel('AoA [deg]'); title('Actual AoA');

%% ========================================================================
%%            PART 3: AERODYNAMIC MAPS (LOOKUP TABLE)
%% ========================================================================
fprintf('\n--- Generating Turbine Performance Maps ---\n');
v_vec   = linspace(10, 30, 50);   % Increased resolution for smoother map
rpm_vec = linspace(1000, 3000, 50); 
[V_GRID, RPM_GRID] = meshgrid(v_vec, rpm_vec);
Torque_Map = zeros(size(V_GRID)); Drag_Map = zeros(size(V_GRID));

for i = 1:size(V_GRID, 1)
    for j = 1:size(V_GRID, 2)
        V_inf = V_GRID(i,j); Omega = RPM_GRID(i,j) * 2 * pi / 60;
        T_sum = 0; D_sum = 0;
        for k = 1:N_span-1
            r_local = (x_span(k) + x_span(k+1))/2;
            dr      = x_span(k+1) - x_span(k);
            c_local = (c_dist(k) + c_dist(k+1))/2;
            
            phi = atan((Omega * r_local) / V_inf);
            beta_local = (beta_phys(k) + beta_phys(k+1))/2;
            alpha_deg = rad2deg(beta_local - phi);
            
            Cl = interp1(Alpha_data, Cl_data, alpha_deg, 'linear', 'extrap');
            Cd = interp1(Alpha_data, Cd_data, alpha_deg, 'linear', 'extrap');
            
            V_rel = sqrt(V_inf^2 + (Omega*r_local)^2);
            L = 0.5 * rho * V_rel^2 * c_local * Cl; 
            D = 0.5 * rho * V_rel^2 * c_local * Cd; 
            
            % TURBINE PHYSICS: Lift drives, Drag resists. Both push backward.
            dF_axial_drag = L * sin(phi) + D * cos(phi);
            dF_driving_torque = L * cos(phi) - D * sin(phi);
            
            T_sum = T_sum + dF_driving_torque * r_local * dr;
            D_sum = D_sum + dF_axial_drag * dr;
        end
        Torque_Map(i,j) = T_sum * B; Drag_Map(i,j) = D_sum * B;
    end
end

%% ============== FIGURE: MAPS ==================
figure('Color','w','Name','Turbine Performance Maps','WindowStyle','docked'); 
subplot(1,2,1); [C,h] = contourf(V_GRID, RPM_GRID, Torque_Map, 40); clabel(C,h); colorbar; colormap(jet); hold on;
% This white line shows your zero-torque freewheeling states!
contour(V_GRID, RPM_GRID, Torque_Map, [0 0], 'w', 'LineWidth', 3);
plot(V_inf_design, RPM_design, 'p', 'MarkerSize', 18, 'MarkerFaceColor',[1 0.84 0], 'MarkerEdgeColor','k');
xlabel('V [m/s]'); ylabel('RPM'); title('TURBINE TORQUE [Nm]');

subplot(1,2,2); [C2,h2] = contourf(V_GRID, RPM_GRID, Drag_Map, 40); clabel(C2,h2); colorbar; colormap(hot); hold on;
plot(V_inf_design, RPM_design, 'p', 'MarkerSize', 18, 'MarkerFaceColor',[1 0.84 0], 'MarkerEdgeColor','k');
xlabel('V [m/s]'); ylabel('RPM'); title('AXIAL DRAG [N]');

%% ========================================================================
%%            PART 4: 3D GEOMETRY PLOT
%% ========================================================================
[sx,sy,sz] = sphere(8); r_marker = c_base/50;
fig = figure('Color','w','Name','Blade Geometry');
axes; hold on; axis equal; grid on
xlabel('X (Radial)'); ylabel('Y (Tangential)'); zlabel('Z (Axial)');
title(['Turbine Geometry (B=' num2str(B) ') - NACA 0012'])
view(35,20); camlight; lighting gouraud
X_air = (X_norm - rotation_position) * c_base; Y_air = Y_norm * c_base;
setappdata(fig,'geom',struct('R_hub',R_hub,'R_duct',R_duct,'hub_h',hub_h,'duct_h',duct_h,'c_blade',c_base,'B',B,'x_span',x_span,'beta',beta_phys,'rotation_position', rotation_position,'X_air',X_air,'Y_air',Y_air,'N_span',N_span,'idx_LE', idx_LE, 'idx_TE', idx_TE,'sx',sx,'sy',sy,'sz',sz,'r_marker',r_marker ));
setappdata(fig,'flags',struct('hub',true,'duct',true,'single',false));
setappdata(fig,'bladeMode',blade_mode);
uicontrol('Style','checkbox','String','Show hub','Value',1,'Units','normalized','Position',[0.02 0.16 0.25 0.05],'Callback',@(s,~)toggleFlag(fig,'hub',s.Value));
uicontrol('Style','checkbox','String','Show duct','Value',1,'Units','normalized','Position',[0.02 0.12 0.25 0.05],'Callback',@(s,~)toggleFlag(fig,'duct',s.Value));
uicontrol('Style','checkbox','String','Single blade only','Value',0,'Units','normalized','Position',[0.02 0.08 0.30 0.05],'Callback',@(s,~)toggleFlag(fig,'single',s.Value));
uicontrol('Style','popupmenu','String',{'Constant Z (Projected)','Constant Chord (Twisted)'},'Value',blade_mode,'Units','normalized','Position',[0.02 0.02 0.2 0.05],'Callback',@(s,~)setBladeMode(fig,s.Value));
redraw(fig);

%% ============== FIGURE 2: TWIST ANGLE ==================
figure('Color','w','Name','Twist Angle'); hold on 
Twist_Plot = rad2deg(beta_phys);
coeffs = polyfit(x_span, Twist_Plot, 1); y_fit_1 = polyval(coeffs, x_span);
coeffs_2 = polyfit(x_span, Twist_Plot, 2); y_fit_2 = polyval(coeffs_2, x_span);
plot(x_span, Twist_Plot, 'LineWidth',2, 'DisplayName', 'Twist Angle');
plot(x_span, y_fit_1, 'DisplayName', 'Linear Fit');
plot(x_span, y_fit_2, 'DisplayName', 'Quadratic Fit');
grid on; legend; xlabel('Span [m]'); ylabel('\beta [deg]'); title('Physical Twist Angle')

%% ============== FIGURE: LE / TE Z vs SPAN ==================
figure('Color','w','Name','LE / TE spanwise');
subplot(1,2,1); grid on; xlabel('Span x [mm]'); ylabel('Z_{LE} [m]'); title('Leading Edge Z (Height)');
subplot(1,2,2); grid on; xlabel('Span x [mm]'); ylabel('Z_{TE} [m]'); title('Trailing Edge Z (Height)');
updateLETEplots(fig, writetable_value);

%% ================= FUNCTIONS =================
function toggleFlag(fig,name,val); F=getappdata(fig,'flags'); F.(name)=logical(val); setappdata(fig,'flags',F); redraw(fig); end
function setBladeMode(fig,val); setappdata(fig,'bladeMode',val); redraw(fig); end
function redraw(fig)
    cla; hold on; G=getappdata(fig,'geom'); F=getappdata(fig,'flags'); mode=getappdata(fig,'bladeMode');
    if F.hub; [th,z]=meshgrid(linspace(0,2*pi,60),linspace(-G.hub_h+G.rotation_position*G.c_blade+1e-3,G.rotation_position*G.c_blade+1e-3,2)); surf(G.R_hub*cos(th),G.R_hub*sin(th),z,'FaceColor',[0.3 0.3 0.3],'EdgeColor','none'); end
    if F.duct; [th,z]=meshgrid(linspace(0,2*pi,80),linspace(-G.duct_h+0.02+G.rotation_position*G.c_blade,0.02+G.rotation_position*G.c_blade,2)); surf(G.R_duct*cos(th),G.R_duct*sin(th),z,'FaceColor',[0.7 0.7 0.7],'EdgeColor','none','FaceAlpha',0.4); end
    beta_root = G.beta(1) + pi; Z_plt_root = G.X_air .* cos(beta_root) + G.Y_air .* sin(beta_root); H_ref = abs(Z_plt_root(G.idx_TE) - Z_plt_root(G.idx_LE));
    if F.single; blades=1; else; blades=1:G.B; end
    for b=blades
        psi=-(b-1)*(2*pi/G.B); Rz=[cos(psi) -sin(psi) 0; sin(psi) cos(psi) 0; 0 0 1];
        for i=1:G.N_span
            beta = G.beta(i) + pi; X_c = G.X_air; Y_t = G.Y_air;
            Y_plt = X_c .* sin(beta) - Y_t .* cos(beta); Z_plt = X_c .* cos(beta) + Y_t .* sin(beta);
            if mode == 1; H_curr = abs(Z_plt(G.idx_TE) - Z_plt(G.idx_LE)); if H_curr > 1e-6; scale = H_ref / H_curr; Y_plt = Y_plt * scale; Z_plt = Z_plt * scale; end; end
            X_plt = zeros(size(Y_plt)); pts = [X_plt; Y_plt; Z_plt]; pts(1,:) = pts(1,:) + G.x_span(i); pts = Rz * pts;
            plot3(pts(1,:), pts(2,:), pts(3,:), 'b');
            if G.idx_LE <= size(pts,2) && G.idx_TE <= size(pts,2)
                LE_pt = pts(:, G.idx_LE); TE_pt = pts(:, G.idx_TE);
                surf(LE_pt(1)+G.r_marker*G.sx, LE_pt(2)+G.r_marker*G.sy, LE_pt(3)+G.r_marker*G.sz, 'FaceColor','r','EdgeColor','none');
                surf(TE_pt(1)+G.r_marker*G.sx, TE_pt(2)+G.r_marker*G.sy, TE_pt(3)+G.r_marker*G.sz, 'FaceColor','g','EdgeColor','none');
            end
        end
    end
    axis equal; grid on; drawnow;
end
function updateLETEplots(fig, want_table)
    G=getappdata(fig,'geom'); mode=getappdata(fig,'bladeMode');
    X_LE=zeros(1,G.N_span); Y_LE=zeros(1,G.N_span); Z_LE=zeros(1,G.N_span); X_TE=zeros(1,G.N_span); Y_TE=zeros(1,G.N_span); Z_TE=zeros(1,G.N_span);
    beta_root = G.beta(1) + pi; Z_plt_root = G.X_air .* cos(beta_root) + G.Y_air .* sin(beta_root); H_ref = abs(Z_plt_root(G.idx_TE) - Z_plt_root(G.idx_LE));
    for i=1:G.N_span
        beta = G.beta(i) + pi; X_c = G.X_air; Y_t = G.Y_air;
        Y_plt = X_c .* sin(beta) - Y_t .* cos(beta); Z_plt = X_c .* cos(beta) + Y_t .* sin(beta);
        if mode == 1; H_curr = abs(Z_plt(G.idx_TE) - Z_plt(G.idx_LE)); if H_curr > 1e-6; scale = H_ref / H_curr; Y_plt = Y_plt * scale; Z_plt = Z_plt * scale; end; end
        X_LE(i) = G.x_span(i); X_TE(i) = G.x_span(i); Y_LE(i) = Y_plt(G.idx_LE); Z_LE(i) = Z_plt(G.idx_LE); Y_TE(i) = Y_plt(G.idx_TE); Z_TE(i) = Z_plt(G.idx_TE);
    end
    if want_table; writetable(table(X_LE', Y_LE', Z_LE', 'VariableNames', {'X_Radial', 'Y_Tangential', 'Z_Axial'}), 'LE_coords.csv'); writetable(table(X_TE', Y_TE', Z_TE', 'VariableNames', {'X_Radial', 'Y_Tangential', 'Z_Axial'}), 'TE_coords.csv'); end
    figure(findobj('Name','LE / TE spanwise')); subplot(1,2,1); cla; plot(G.x_span*1e3, Z_LE, 'LineWidth',2); grid on; xlabel('Span [mm]'); title('LE Z (Height)'); subplot(1,2,2); cla; plot(G.x_span*1e3, Z_TE, 'LineWidth',2); grid on; xlabel('Span [mm]'); title('TE Z (Height)');
end