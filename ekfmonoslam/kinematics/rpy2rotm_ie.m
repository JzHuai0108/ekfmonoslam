function ll_R_b = rpy2rotm_ie(rpy)
% Compute the rotation matrix from the body frame to the local level frame
% in the Inertial Explorer convention.
% The body frame is in general front forward up.
% The local level frame is east north up.
% ll_R_b = R3(-yaw) * R1(-pitch) * R2(-roll). Note our R1(a) = IE R1(-a).
% see eq 6 page 233 in Waypoint Software 8.90 User Manual v10 Aug 2020.
r = rpy(1);
p = rpy(2);
y = rpy(3);
cr = cos(r);
sr = sin(r);
cp = cos(p);
sp = sin(p);
cy = cos(y);
sy = sin(y);
ll_R_b = [cy * cr - sy * sp * sr, -sy * cp, cy * sr + sy * sp * cr;
    sy * cr + cy * sp * sr,  cy * cp, sy * sr - cy * sp * cr;
    -cp * sr, sp,  cp * cr];
% the above is the same as 
% ll_R_b = R3(-y) * R1(-p) * R2(-r);
end