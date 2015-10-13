%otyp==0->cont, otyp==1->disc
function [A N]=sys_bd_dcm(llh, vel_b, Cbn, gyro, otyp, dt)
% position   (1-3)
% velocity   (4-6)
% attitude   (7-9) 

%system disturbance coefs
N=zeros(9,6);
N(7:9,4:6)=-Cbn; %attitude
N(4:6,1:3)=eye(3); %velocity
N(4:6,4:6)=skew(vel_b); %velocity

%system matrix
A=zeros(9);
[Rn, Re, g, sL, cL, WIE_E]=geoparam_v000(llh);
tL=sL/cL;
Rn_h=Rn+llh(3);
Re_h=Re+llh(3);
R=diag([1/Rn_h,1/Re_h,-1]);
O=[0 1/Re_h 0;-1/Rn_h 0 0;0 -tL/Re_h 0];
vel_n=Cbn*vel_b;
wen_n=[vel_n(2)/Re_h; -vel_n(1)/Rn_h; -vel_n(2)*tL/Re_h];
wie_n=[WIE_E*cL; 0; -WIE_E*sL];
win_n=wen_n+wie_n;

%Position
A(1:3,4:6)=R*Cbn;
A(1:3,7:9)=R*skew(vel_n);
A(1,3)=-vel_n(1)/Rn_h^2;
A(2,1)=vel_n(2)*tL/Re_h/cL;
A(2,3)=-vel_n(2)/Re_h^2/cL;

%Velocity
A(4:6,4:6)=-skew(Cbn'*wie_n+gyro);
A(4:6,7:9)=Cbn'*skew([0;0;-g])+skew(vel_b)*Cbn'*skew(wie_n);
A(4:6,1)=skew(vel_b)*Cbn'*[-WIE_E*sL;0;-WIE_E*cL];

%Attitude
A(7,1)=-WIE_E*sL;
A(7,3)=-vel_n(2)/Re_h^2;
A(8,3)=vel_n(1)/Rn_h^2;
A(9,1)=-(WIE_E*cL+vel_n(2)/Re_h/cL/cL);
A(9,3)=vel_n(2)*tL/Re_h^2;
A(7:9,4:6)=O*Cbn;
A(7:9,7:9)=-skew(win_n)+O*skew(vel_n);

%discretize A;
if (otyp==0)
    A=A;
elseif(otyp==1)
    A=expm(A*dt);
end


