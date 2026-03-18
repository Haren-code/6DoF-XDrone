function q = quatNorm(q_original)
%QUATNORM Normalize quaternion (scalar-first convention)
% Explicit multiplication is faster than exponentiation
    norm_q2 = q_original(1)*q_original(1) + ...
              q_original(2)*q_original(2) + ...
              q_original(3)*q_original(3) + ...
              q_original(4)*q_original(4); 

    if norm_q2 > 1e-24
        % Multiply by inverse square root
        q = q_original * (1 / sqrt(norm_q2)); 
    else
        q = [1;0;0;0]; % identity quaternion
        warning('quatNorm:ZeroQuaternion','Quaternion norm nearly zero. Reset to identity.');
    end
    
end