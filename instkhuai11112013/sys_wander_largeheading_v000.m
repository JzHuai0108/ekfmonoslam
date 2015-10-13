%System model for large initial uncertainty (wander-azimuth mechanization,
%attitude=dcm, pos=llh). The large wander angle is modelled as [sin(alpha), cos(alpha)]
%The error definitions are based on phi-model
%Depending on input arguments this function may return (A,N), (Ad,N),
%(Ac,Qc),(Ad,Qd)
%mode==1->fine alignment, otherwise coarse alignment

%(For an example usage of this see: example_largeheading)

function [RetA RetB]=sys_wander_largeheading(llh, vel_n, Cbn, wander, acc, dt, mode, Aimu, Qimu, Cimu, Rimu)
% position   (1-3) (not it llh but in northing, easting, height)
% velocity   (4-6)
% attitude   (7-9)
% wander     (10-11) ([del_sa, del_ca])

%system disturbance coefs
N=zeros(11,6);
N(7:9,4:6)=-Cbn; %attitude
N(4:6,1:3)=Cbn; %velocity


%system matrix
A=zeros(11);

[Rn, Re, g, sL, cL, WIE_E]=geoparam_v000(llh);  %do not need to be exact by any means.
Cgn=[wander(2) wander(1) 0;-wander(1) wander(2) 0; 0 0 1];
Rn_h=Rn+llh(3);
Re_h=Re+llh(3);
R=sqrt(Rn_h*Re_h);
acc_n=Cbn*acc;

mx_a=[0 WIE_E*cL;-WIE_E*sL 0;0 0];
mx_b=[wander(2)*(-WIE_E)*sL/Rn_h 0;wander(1)*(-WIE_E)*sL/Rn_h 0; (-WIE_E)*cL/Rn_h 0];
mx_c=[0 1/R 0;-1/R 0 0;0 0 0];
mx_d=[-vel_n(2)/R^2; vel_n(1)/R^2;0];


mx_e=Cgn*[0 1/(Re_h) 0; -1/(Rn_h) 0 0;0 0 0]*Cgn'*vel_n+2*Cgn*[WIE_E*cL; 0; -WIE_E*sL]; %2wie_n+wen_n
mx_f=Cgn*[0 1/(Re_h) 0; -1/(Rn_h) 0 0;0 0 0]*Cgn'*vel_n+Cgn*[WIE_E*cL; 0; -WIE_E*sL]; %wie_n+wen_n

%%position errors
A(1:2,4:5)=Cgn(1:2,1:2)';
A(1:2,10:11)=[-vel_n(2) vel_n(1);vel_n(1) vel_n(2)];
A(3,6)=-1;

%%Velocity errors
A(4:6,1:2)=2*skew(vel_n)*mx_b;
A(4:6,3)=skew(vel_n)*mx_d+[0;0;-2*g/R];
A(4:6,4:6)=skew(vel_n)*mx_c-skew(mx_e);
A(4:6,7:9)=skew(acc_n);
A(4:6,10:11)=2*skew(vel_n)*mx_a;


%%Attitude errors
A(7:9,1:2)=mx_b;
A(7:9,3)=mx_d;
A(7:9,4:6)=mx_c;
A(7:9,7:9)=-skew(mx_f);
A(7:9,10:11)=mx_a;

%%Wander Errors
A(10:11,10:11)=0;


if (mode==1)    %convert the model for fine alignment
    mx_a=zeros(10,10);
    mx_a(:,1:9)=A(1:10,1:9);
    mx_a(:,10)=wander(2)*A(1:10,10)-wander(1)*A(1:10,11);
    A=mx_a;
    N=N(1:10,:);
end
    


%%%%%%%%
%%Combine navigation models with imu models and return
if ~(isempty(Aimu) || isempty(Qimu) || isempty(Rimu) || isempty(Cimu))
    nst_imu=size(Aimu,1);
    nst_nav=size(A,1);
    Ac=[A N*Cimu;zeros(nst_imu,nst_nav) Aimu];
    Qc=[N*Rimu*N' zeros(nst_nav,nst_imu);zeros(nst_imu,nst_nav) Qimu];
    
    RetA=Ac;
    RetB=Qc;
    
    if (dt>0)
        %Discretize the model
        nst=nst_nav+nst_imu;
        mx_a = dt*[-Ac,Qc;zeros(nst),Ac'];
        mx_b = expm(mx_a);
        STM = mx_b(nst+1:2*nst,nst+1:2*nst)';
        Qd = STM*mx_b(1:nst,nst+1:2*nst);
        RetA=STM;
        RetB=Qd;
    end
else
    RetA=A;
    RetB=N;
    if (dt>0)
        RetA=expm(A*dt);
    end
end
