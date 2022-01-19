function rpy = rotm2rpy_ie(ll_R_b)
% Compute the [roll; pitch; yaw] in radians given rotation matrix from 
% the body frame to the local level frame in the Inertial Explorer convention.
% The body frame is in general front forward up.
% The local level frame is east north up.
% ll_R_b = R3(-yaw) * R1(-pitch) * R2(-roll). Note our R1(a) = IE R1(-a).
% see eq 6 page 233 in Waypoint Software 8.90 User Manual v10 Aug 2020.
r = atan2(-ll_R_b(3, 1), ll_R_b(3, 3));
p = asin(ll_R_b(3, 2));
y = atan2(-ll_R_b(1, 2), ll_R_b(2, 2));
rpy = [r; p; y];
end
