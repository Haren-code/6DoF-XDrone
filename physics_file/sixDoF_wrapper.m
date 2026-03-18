function dx = sixDoF_wrapper(t, x, sim)
    
    % Unpack vector x containing all the states
    v_b = x(1:3);        % First 3 states are velocities
    omega_b = x(4:6);    % Next 3 states are angular velocities
    q = x(7:10);         % Next 4 states are quaternion parameters
    %pos_i = x(11:13);   % Last 3 are the position vector - bot not used here
    
    % Normalize quaternion
    q = quatNorm(q);
    
    % Compute rotation matrices for rotating body frame to inertial frame
    R_bi = quatRotMat(q);    % quaterinon to R_bi = body -> inertial
    R_ib = R_bi';            % R_ib = inertial -> body
      
    % Dynamics (calculates body forces and moments)
    [F_b, M_b] = dynamics(t, v_b, omega_b, R_ib, sim); 
    
    % Kinematics (calculates derivatives of the states)
    [v_dot_b, omega_dot_b, q_dot, pos_dot_i] = kinematics(v_b, omega_b, F_b, M_b, q, R_bi, sim);
    
    % Pack state derivatives vector for integration (main.m has the integration function)
    dx = [v_dot_b; omega_dot_b; q_dot; pos_dot_i];

end