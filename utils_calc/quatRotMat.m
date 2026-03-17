function R = quatRotMat(q)
% QUATROTMAT Converts quaternions to rotational matrix
% Assumes q is already normalized, scalar-first [w x y z]
    w = q(1); x = q(2); y = q(3); z = q(4);
    
    x2 = x*x; y2 = y*y; z2 = z*z;
    wx = w*x; wy = w*y; wz = w*z;
    xy = x*y; xz = x*z; yz = y*z;
    
    R = [1-2*(y2+z2),  2*(xy-wz),   2*(xz+wy);
           2*(xy+wz),  1-2*(x2+z2), 2*(yz-wx);
           2*(xz-wy),  2*(yz+wx),   1-2*(x2+y2)];
end