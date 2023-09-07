function ll_R_b = rpy2rotm_ie(rpy)
% Compute the rotation matrix from the body frame to the local level frame
% in the Inertial Explorer convention.
% The body frame is in general right forward up.
% The local level frame is east north up.
% ll_R_b = R3(-yaw) * R1(-pitch) * R2(-roll). Note our R1(a) = IE R1(-a).
% see page 50 bottom in waypoint software 6.0 user manual November 2014
% and see eq 6 page 233 in Waypoint Software 8.90 User Manual v10 Aug 2020.
% Note bynav uses the same notation as described here
% https://www.bynav.com/cn/resource/bywork/healthy-work/654.html
% About heading/azimuth and yaw, according to Waypoint Software 8.90 User
% Manual v10 Aug 2020, "Positive heading rotation is clockwise from North",
% and "yaw = –heading". As the manual plots azimuth and heading together,
% I think, azimuth is equivalent to heading. This is also confirmed here
% https://aviation.stackexchange.com/questions/24899/what-is-the-difference-between-azimuth-and-heading

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