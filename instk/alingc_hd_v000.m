%Let Cpn_est=(I-S(e))Cpn. Then e=[Egyro Eplat]*[gyro error;Platform error]
%platform err def:Cbp_est=(I-S(p))*Cbp
%gyro err def: gyro_meas=gyro+gyro_error
function [Cnp, Egyro, Eplat]=alingc_hd(gyro, Lat, Cbp)

WIE_E=7292115e-11;
nrot=WIE_E*cos(Lat);

%normalize gyro;
gyro=gyro/norm(gyro)*WIE_E;

%Compute the earth rotation in platform frame
if isempty(Cbp)
    Cbp=eye(3);
end
gyro_p=Cbp*gyro;

%Transformation
Cnp=zeros(3);
Cnp(1,1)=gyro_p(1)/nrot;
Cnp(2,2)=gyro_p(1)/nrot;
Cnp(2,1)=gyro_p(2)/nrot;
Cnp(1,2)=-Cnp(2,1);
Cnp(3,3)=1;

%Error matrices
Egyro=[-gyro_p(2)/nrot^2 gyro_p(1)/nrot^2 0]*Cbp;
Eplat=[-gyro_p(2)/nrot^2 gyro_p(1)/nrot^2 0]*skew(gyro_p);
