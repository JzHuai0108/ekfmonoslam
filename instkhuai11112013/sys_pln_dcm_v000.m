%vtyp==0->n frame, vtyp==1->b frame
%otyp==0->cont, otyp==1->disc
function [A N]=sys_pln_dcm(pos_n, vel_x, Cbn, acc, gyro, g, vtyp, otyp, dt)
% position   (1-3)
% velocity   (4-6)
% attitude   (7-9) 

A=zeros(9);
N=zeros(9,6);
gravity=[0;0;-g];

%%%IMU error inputs
%attitude
N(7:9,4:6)=-Cbn;
%velocity
if (vtyp==0)
    N(4:6,1:3)=Cbn;
elseif (vtyp==1)
    N(4:6,1:3)=eye(3);
    N(4:6,4:6)=skew(vel_x);
end
%position (no direct imu effect)


%%%Continuous time system parameters
F=zeros(9);
%atttitude (no dynamic)
%Velocity
if (vtyp==0)
    F(4:6,7:9)=skew(Cbn*acc);
elseif (vtyp==1)
    F(4:6,4:6)=-skew(gyro);
    F(4:6,7:9)=Cbn'*skew(gravity);      %%%it hurts
end
%position
if (vtyp==0)
    F(1:3,4:6)=eye(3);
elseif (vtyp==1)
    F(1:3,4:6)=Cbn;
    F(1:3,7:9)=skew(Cbn*vel_x);
end

%discretize F;
if (otyp==0)
    A=F;
elseif(otyp==1)
    A=eye(9)+F*dt+F*F*dt*dt/2;
end