%Given measured acc and any Azimuth estimate, computes the nav to body frame rotation
%Cpn=euler2dcm([0;0;AZ]).
%Let Cbp_est=(I-S(e1))Cbp. Then e1=E*(acc error).
%Let Cbn_est=(I-S(e2))Cbn. Then e2=e1+[0;0; (AZ error)].

function [Cnb E]=alingc_pl(acc, AZ)
g=norm(acc);
acc=acc/g; %normalize
den=sqrt(1-(acc(1))^2);    %Note:assumes pitch is not 90

%Platform to body transformation
Cpb=zeros(3,3);
Cpb(:,3)=-acc;

Cpb(2,2)=-acc(3)/den;
Cpb(3,2)=acc(2)/den;

Cpb(1,1)=den;
Cpb(2,1)=-acc(1)*acc(2)/den;
Cpb(3,1)=-acc(1)*acc(3)/den;

%Error matrix
E=zeros(3);
E(1,2)=-acc(3)/den/g;
E(1,3)=acc(2)/den/g;

E(2,1)=-den/g;
E(2,2)=acc(1)*acc(2)/den/g;
E(2,3)=acc(1)*acc(3)/den/g;

E(3,2)=acc(1)*acc(3)/den^2/g;
E(3,3)=-acc(1)*acc(2)/den^2/g;

%platform to NED
if (isempty(AZ))
    AZ=0;
end
Cpn=euler2dcm_v000([0;0;AZ]);
    
%NED to Body
Cnb=Cpb*Cpn';
E=Cpn*E;

