%Standard yet most useless method.
function Cnb=alingc_std(acc,gyro,Llh, mtd)

[Rn, Re, g, sL, cL, WIE_E]=geoparam_v000(Llh);

if (mtd==0)
    %normalize acc ang gyro;
    acc=acc/norm(acc);
    gyro=gyro/norm(gyro);
    Cnb=[cross(acc,gyro) gyro acc]*[0 -1/cL 0;1/cL 0 0;-sL/cL 0 -1];
elseif (mtd==1) %Rogers, Titterton etc
    Cnb=[cross(acc,gyro) gyro acc]*[0 -1/cL/WIE_E/g 0;1/cL/WIE_E 0 0;-sL/cL/g 0 -1/g];
end
