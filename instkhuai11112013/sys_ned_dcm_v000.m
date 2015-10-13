%otyp==0->cont, otyp==1->disc
function [A N]=sys_ned_dcm(llh, vel_n, Cbn, acc, otyp, dt)
% position   (1-3)
% velocity   (4-6)
% attitude   (7-9) 

%system disturbance coefs
N=zeros(9,6);
N(7:9,4:6)=-Cbn; %attitude
N(4:6,1:3)=Cbn; %velocity

%system matrix
A=zeros(9);
[Rn, Re, g, sL, cL, WIE_E]=geoparam_v000(llh);
tL=sL/cL;
Rn_h=Rn+llh(3);
Re_h=Re+llh(3);
acc_n=Cbn*acc;

A(1,3)=-vel_n(1)/(Rn_h)^2; 
A(1,4)=1/(Rn_h);

A(2,1)=vel_n(2)*tL/(Re_h)/cL;
A(2,3)=-vel_n(2)/(Re_h)/(Re_h)/cL;
A(2,5)=1/(Re_h)/cL;

A(3,6)=-1;

A(4,1)=-(2*WIE_E*cL*vel_n(2)+vel_n(2)*vel_n(2)/(Re_h)/cL/cL);
A(4,3)=(-vel_n(1)*vel_n(3)/(Rn_h)/(Rn_h)+vel_n(2)*vel_n(2)*tL/(Re_h)/(Re_h));
A(4,4)=vel_n(3)/(Rn_h);
A(4,5)=-(2*WIE_E*sL+2*vel_n(2)*tL/(Re_h));
A(4,6)=vel_n(1)/(Rn_h);
A(4,8)=-acc_n(3);
A(4,9)=acc_n(2);

A(5,1)=(2*WIE_E*cL*vel_n(1)-2*WIE_E*sL*vel_n(3)+vel_n(2)*vel_n(1)/(Re_h)/cL/cL);
A(5,3)=-(vel_n(2)*vel_n(3)+vel_n(2)*vel_n(1)*tL)/(Re_h)/(Re_h);
A(5,4)=(2*WIE_E*sL+vel_n(2)*tL/(Re_h));
A(5,5)=(vel_n(3)+vel_n(1)*tL)/(Re_h);
A(5,6)=(2*WIE_E*cL+vel_n(2)/(Re_h));
A(5,7)=acc_n(3);
A(5,9)=-acc_n(1);

A(6,1)=2*WIE_E*sL*vel_n(2);
A(6,3)=(vel_n(1)*vel_n(1)/(Rn_h)/(Rn_h)+vel_n(2)*vel_n(2)/(Re_h)/(Re_h));
A(6,4)=-2*vel_n(1)/(Rn_h);
A(6,5)=-2*(WIE_E*cL+vel_n(2)/(Re_h));
A(6,7)=-acc_n(2);
A(6,8)=acc_n(1);


A(7,1)=-WIE_E*sL;
A(7,3)=-vel_n(2)/(Re_h)/(Re_h);
A(7,5)=1/(Re_h);
A(7,8)=-(WIE_E*sL+vel_n(2)*tL/(Re_h));
A(7,9)=vel_n(1)/(Rn_h);

A(8,3)=vel_n(1)/(Rn_h)/(Rn_h);
A(8,4)=-1/(Rn_h);
A(8,7)=(WIE_E*sL+vel_n(2)*tL/(Re_h));
A(8,9)=(WIE_E*cL+vel_n(2)/(Re_h));

A(9,1)=-(WIE_E*cL+vel_n(2)/(Re_h)/cL/cL);
A(9,3)=vel_n(2)*tL/(Re_h)/(Re_h);
A(9,5)=-tL/(Re_h);
A(9,7)=-vel_n(1)/(Rn_h);
A(9,8)=-(WIE_E*cL+vel_n(2)/(Re_h));


%discretize A;
if (otyp==0)
    A=A;
elseif(otyp==1)
    A=expm(A*dt);
end


