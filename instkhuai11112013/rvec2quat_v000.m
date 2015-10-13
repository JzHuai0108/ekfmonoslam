function quat=rvec2quat(rvec)

%normalize
rot_ang=sqrt(rvec(1)^2+rvec(2)^2+rvec(3)^2);    %assume always positive
if (rot_ang==0)
    quat=[1;0;0;0];
else
    cR=cos(rot_ang/2);
    sR=sin(rot_ang/2);
    rx=rvec(1)/rot_ang;
    ry=rvec(2)/rot_ang;
    rz=rvec(3)/rot_ang;

    quat=[cR;rx*sR;ry*sR;rz*sR];
end
