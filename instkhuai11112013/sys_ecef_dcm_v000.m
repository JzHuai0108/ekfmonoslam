%otyp==0->cont, otyp==1->disc
function [A N]=sys_ecef_dcm(Cbe, acc, otyp, dt)
% position   (1-3)
% velocity   (4-6)
% attitude   (7-9) 

WIE_E=7292115e-11;  %Earth rotation rate

%system disturbance coefs
N=zeros(9,6);
N(7:9,4:6)=-Cbe; %attitude
N(4:6,1:3)=Cbe; %velocity

%system matrix
A=zeros(9);
acc_e=Cbe*acc;

%Position
A(1,4)=1;
A(2,5)=1;
A(3,6)=1;

%Velocity
A(4:6,4:6)=-2*skew([0;0;WIE_E]);
A(4:6,7:9)=skew(acc_e);

%Attitude
A(7:9,7:9)=-skew([0;0;WIE_E]);

%discretize A;
if (otyp==0)
    A=A;
elseif(otyp==1)
    A=expm(A*dt);
end


