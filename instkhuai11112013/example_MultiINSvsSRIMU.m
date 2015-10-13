clear;

%main nav states
nst_n=2;
An=diag([1 1]);
N=eye(2);
sPn=eye(2)*0.01;
xn=sPn*randn(nst_n,1);

%base sensor model
[At, Qt, Rt]=dc2dc_v000(0.99998,(0.5e-5)^2,0.02^2,1/100,2,1);

%Sensors
A1=At;
B1=sqrt(Qt);
C1=[1];
M1=[cos(pi/11) sin(pi/11)];
sP1=sqrt(B1^2/(1-A1^2));
x1=sP1*randn(1,1);
sR1=sqrt(Rt);

A2=At;
B2=sqrt(Qt*10);
C2=[1];
M2=[cos(2*pi/11) sin(2*pi/11)];
sP2=sqrt(B2^2/(1-A2^2));
x2=sP2*randn(1,1);
sR2=sqrt(Rt/10);

A3=At;
B3=sqrt(Qt*100);
C3=[1];
M3=[cos(3*pi/11) sin(3*pi/11)];
sP3=sqrt(B3^2/(1-A3^2));
x3=sP3*randn(1,1);
sR3=sqrt(Rt/100);

A4=At;
B4=sqrt(Qt/10);
C4=[1];
M4=[cos(4*pi/11) sin(4*pi/11)];
sP4=sqrt(B4^2/(1-A4^2));
x4=sP4*randn(1,1);
sR4=sqrt(Rt*10);

A5=At;
B5=sqrt(Qt/100);
C5=[1];
M5=[cos(5*pi/11) sin(5*pi/11)];
sP5=sqrt(B5^2/(1-A5^2));
x5=sP5*randn(1,1);
sR5=sqrt(Rt*100);

%%% Main model
%imu observation structure
Rimu=sR1*sR1';
Rimu=diagmat_v000(sR2*sR2',Rimu,1);
Rimu=diagmat_v000(sR3*sR3',Rimu,1);
Rimu=diagmat_v000(sR4*sR4',Rimu,1);
Rimu=diagmat_v000(sR5*sR5',Rimu,1);

Cimu=C1;
Cimu=diagmat_v000(C2,Cimu,1);
Cimu=diagmat_v000(C3,Cimu,1);
Cimu=diagmat_v000(C4,Cimu,1);
Cimu=diagmat_v000(C5,Cimu,1);

M=[M1;M2;M3;M4;M5];

%single observation model and input
[T Mls]=cp_Tparam_v000(M,Rimu);
Cimu_obs=[zeros(size(T,1), nst_n) T*Cimu];
Rimu_obs=T*Rimu*T';
Cimu_inp=Mls*Cimu;
Rimu_inp=Mls*Rimu*Mls';

%single system model
A=An;
A=diagmat_v000(A1,A,1);
A=diagmat_v000(A2,A,1);
A=diagmat_v000(A3,A,1);
A=diagmat_v000(A4,A,1);
A=diagmat_v000(A5,A,1);
A(1:nst_n,(nst_n+1):end)=-N*Cimu_inp;

Q=N*Rimu_inp*N';
Q=diagmat_v000(B1*B1',Q,1);
Q=diagmat_v000(B2*B2',Q,1);
Q=diagmat_v000(B3*B3',Q,1);
Q=diagmat_v000(B4*B4',Q,1);
Q=diagmat_v000(B5*B5',Q,1);

P=sPn*sPn';
P=diagmat_v000(sP1*sP1',P,1);
P=diagmat_v000(sP2*sP2',P,1);
P=diagmat_v000(sP3*sP3',P,1);
P=diagmat_v000(sP4*sP4',P,1);
P=diagmat_v000(sP5*sP5',P,1);

xest=zeros(size(A,1),1);


%%%%% 2-sub system model
%observation structure
Ra=sR1*sR1';
Ra=diagmat_v000(sR2*sR2',Ra,1);
Ra=diagmat_v000(sR3*sR3',Ra,1);
Rb=sR4*sR4';
Rb=diagmat_v000(sR5*sR5',Rb,1);

Cimua=C1;
Cimua=diagmat_v000(C2,Cimua,1);
Cimua=diagmat_v000(C3,Cimua,1);
Cimub=C4;
Cimub=diagmat_v000(C5,Cimub,1);

Ma=[M1;M2;M3];
Mb=[M4;M5];

%observation models and inputs
[Ta Mlsa]=cp_Tparam_v000(Ma,Ra);
Mlsb=lsmat_v001(Mb,Rb,0);

Cimu_obsa=[zeros(size(Ta,1), 2*nst_n) Ta*Cimua zeros(size(Ta,1), size(Cimub,2))];
Rimu_obsa=Ta*Ra*Ta';

Cimu_inpa=Mlsa*Cimua;
Rimu_inpa=Mlsa*Ra*Mlsa';
Cimu_inpb=Mlsb*Cimub;
Rimu_inpb=Mlsb*Rb*Mlsb';

%system models
Aab=An;
Aab=diagmat_v000(An,Aab,1);
Aab=diagmat_v000(A1,Aab,1);
Aab=diagmat_v000(A2,Aab,1);
Aab=diagmat_v000(A3,Aab,1);
Aab=diagmat_v000(A4,Aab,1);
Aab=diagmat_v000(A5,Aab,1);
Aab(1:nst_n,(2*nst_n+(1:size(Cimu_inpa,2))))=-N*Cimu_inpa;
Aab(nst_n+1:2*nst_n,((2*nst_n+size(Cimu_inpa,2))+(1:size(Cimu_inpb,2))))=-N*Cimu_inpb;

Qab=N*Rimu_inpa*N';
Qab=diagmat_v000(N*Rimu_inpb*N',Qab,1);
Qab=diagmat_v000(B1*B1',Qab,1);
Qab=diagmat_v000(B2*B2',Qab,1);
Qab=diagmat_v000(B3*B3',Qab,1);
Qab=diagmat_v000(B4*B4',Qab,1);
Qab=diagmat_v000(B5*B5',Qab,1);

Pab=sPn*sPn';
Pab=diagmat_v000(Pab,Pab,1);
Pab=diagmat_v000(sP1*sP1',Pab,1);
Pab=diagmat_v000(sP2*sP2',Pab,1);
Pab=diagmat_v000(sP3*sP3',Pab,1);
Pab=diagmat_v000(sP4*sP4',Pab,1);
Pab=diagmat_v000(sP5*sP5',Pab,1);
Pab(1:nst_n,nst_n+1:2*nst_n)=sPn*sPn';
Pab(nst_n+1:2*nst_n,1:nst_n)=sPn*sPn';

xest_ab=zeros(size(Aab,1),1);

%cross_correlation corrected system
mx_d=zeros(size(Ra,1), size(Rimu,2));
mx_d(:,1:size(Ra,1))=Ra;
mx_c=zeros(size(Rb,1), size(Rimu,2));
mx_c(:,(size(Ra,1)+1):end)=Rb;
mx_a=-N*Mlsa*mx_d*T';
mx_b=-N*Mlsb*mx_c*T';
S=zeros(size(Aab,1),size(T,1));
S(1:nst_n,:)=mx_a;
S((nst_n+1):2*nst_n,:)=mx_b;

Cimu_obs_cr=[zeros(size(T,1), 2*nst_n) T*Cimu];

mx_a=inv(Rimu_obs);
Aab_cr=Aab-S*mx_a*Cimu_obs_cr;
Pab_cr=Pab;
Qab_cr=Qab-S*mx_a*S';
N_cr=S*mx_a;

xest_ab_cr=zeros(size(Aab,1),1);


%system under state equivalence
Aab_eq=Aab;
Pab_eq=Pab;
Qab_eq=Qab;
Cimu_obs_eq=Cimu_obsa;
Cst_eq=zeros(nst_n, size(Aab_eq,1));
Cst_eq(1:nst_n,1:nst_n)=eye(nst_n);
Cst_eq(1:nst_n,(nst_n+1):2*nst_n)=-eye(nst_n);
Rimu_obs_eq=Rimu_obsa;
xest_ab_eq=zeros(size(Aab,1),1);

%system under state equivalence 1
Aab_eq1=Aab;
Pab_eq1=Pab;
Qab_eq1=Qab;
Cst_eq1=zeros(nst_n, size(Aab_eq,1));
Cst_eq1(1:nst_n,1:nst_n)=eye(nst_n);
Cst_eq1(1:nst_n,(nst_n+1):2*nst_n)=-eye(nst_n);
xest_ab_eq1=zeros(size(Aab,1),1);


%system under state equivalence 2
Aab_eq2=Aab;
Pab_eq2=Pab;
Qab_eq2=Qab;

mx_a=[zeros(size(Ta,1), 2*nst_n) Ta*Cimua zeros(size(Ta,1), size(Cimub,2))]; %observation of redundant sensor in set a
mx_c=[zeros(size(Cimu_inpa,1), 2*nst_n) Cimu_inpa -Cimu_inpb]; %observation of seta-setb
Cimu_obs_eq2=[mx_a;mx_c];
mx_b=[Ta zeros(size(Ta,1),size(Mlsb,2));[Mlsa -Mlsb]];
Rimu_obs_eq2=mx_b*Rimu*mx_b';

xest_ab_eq2=zeros(size(Aab,1),1);


%system under state equivalence 3
Aab_eq3=Aab;
Pab_eq3=Pab;
Qab_eq3=Qab;

mx_a=[zeros(size(Ta,1), 2*nst_n) Ta*Cimua zeros(size(Ta,1), size(Cimub,2))]; %observation of redundant sensor in set a
mx_c=[zeros(size(Cimu_inpa,1), 2*nst_n) Cimu_inpa -Cimu_inpb]; %observation of seta-setb
Cimu_obs_eq3=[mx_c];

mx_b=[Ta zeros(size(Ta,1),size(Mlsb,2));[Mlsa -Mlsb]];
Rimu_obs_eq3=[Mlsa -Mlsb]*Rimu*[Mlsa -Mlsb]';

xest_ab_eq3=zeros(size(Aab,1),1);

%%%%%%%%%%%%%%%%%%
nit=300;
nstab=size(Aab,1);
nst=size(A,1);
%debug variables
debctr=0;
debx=zeros(nst,nit);        %sol 1
debxval=zeros(nst,nit);     %true vals
debp=zeros(nst,nit);
debxab=zeros(nstab,nit);    %not used
debpab=zeros(nstab,nit);
debxabcr=zeros(nstab,nit);  %sol 5a
debpabcr=zeros(nstab,nit);
debxeq=zeros(nstab,nit);    %sol 3
debpeq=zeros(nstab,nit); 
debxeq1=zeros(nstab,nit);   %sol 2
debpeq1=zeros(nstab,nit);
debxeq2=zeros(nstab,nit);   %sol 5b
debpeq2=zeros(nstab,nit);
debxeq3=zeros(nstab,nit);   %sol 4
debpeq3=zeros(nstab,nit);

debctr=debctr+1;
debx(:,debctr)=xest;
debxval(:,debctr)=[xn;x1;x2;x3;x4;x5];
debp(:,debctr)=diag(P).^0.5;
debxab(:,debctr)=xest_ab;
debpab(:,debctr)=diag(Pab).^0.5;
debxabcr(:,debctr)=xest_ab_cr;
debpabcr(:,debctr)=diag(Pab_cr).^0.5;
debxeq(:,debctr)=xest_ab_eq;
debpeq(:,debctr)=diag(Pab_eq).^0.5;
debxeq1(:,debctr)=xest_ab_eq1;
debpeq1(:,debctr)=diag(Pab_eq1).^0.5;
debxeq2(:,debctr)=xest_ab_eq2;
debpeq2(:,debctr)=diag(Pab_eq2).^0.5;
debxeq3(:,debctr)=xest_ab_eq3;
debpeq3(:,debctr)=diag(Pab_eq3).^0.5;

%Start main loop
randn('state',243243);
flag=1;
for in=1:nit
    %observations
    u=randn(2,1);
    y1=M1*u+C1*x1+sR1*randn(size(sR1,1),1);
    y2=M2*u+C2*x2+sR2*randn(size(sR2,1),1);
    y3=M3*u+C3*x3+sR3*randn(size(sR3,1),1);
    y4=M4*u+C4*x4+sR4*randn(size(sR4,1),1);
    y5=M5*u+C5*x5+sR5*randn(size(sR5,1),1);
   
    y=[y1;y2;y3;y4;y5];
    imu_inp=Mls*y;
    imu_obs=T*y;
    
    y=[y1;y2;y3];
    imu_inpa=Mlsa*y;
    imu_obsa=Ta*y;
    
    y=[y4;y5];
    imu_inpb=Mlsb*y;
    
    imu_obs_eq2=imu_inpa-imu_inpb;
    
    
    %Kalman optimal
    K=P*Cimu_obs'*inv(Cimu_obs*P*Cimu_obs'+Rimu_obs);
    P=P-K*Cimu_obs*P;
    dz=K*(imu_obs-Cimu_obs*xest);
    xest=xest+dz;
      
    %Kalman ab
    Kab=Pab*Cimu_obsa'*inv(Cimu_obsa*Pab*Cimu_obsa'+Rimu_obsa);
    Pab=Pab-Kab*Cimu_obsa*Pab;
    dz=Kab*(imu_obsa-Cimu_obsa*xest_ab);
    xest_ab=xest_ab+dz;
    
    %kalman ab cross_correlation corrected
    Kab_cr=Pab_cr*Cimu_obs_cr'*inv(Cimu_obs_cr*Pab_cr*Cimu_obs_cr'+Rimu_obs);
    Pab_cr=Pab_cr-Kab_cr*Cimu_obs_cr*Pab_cr;
    dz=Kab_cr*(imu_obs-Cimu_obs_cr*xest_ab_cr);
    xest_ab_cr=xest_ab_cr+dz;
    
    %kalman equivalence
    Kab_eq=Pab_eq*Cimu_obs_eq'*inv(Cimu_obs_eq*Pab_eq*Cimu_obs_eq'+Rimu_obs_eq);
    Pab_eq=Pab_eq-Kab_eq*Cimu_obs_eq*Pab_eq;
    dz=Kab_eq*(imu_obsa-Cimu_obs_eq*xest_ab_eq);
    xest_ab_eq=xest_ab_eq+dz;
    
    %Kalman equivalence2
    Kab_eq2=Pab_eq2*Cimu_obs_eq2'*inv(Cimu_obs_eq2*Pab_eq2*Cimu_obs_eq2'+Rimu_obs_eq2);
    Pab_eq2=Pab_eq2-Kab_eq2*Cimu_obs_eq2*Pab_eq2;
    dz=Kab_eq2*([imu_obsa;imu_obs_eq2]-Cimu_obs_eq2*xest_ab_eq2);
    xest_ab_eq2=xest_ab_eq2+dz;
    
    
    %Kalman equivalence3
    Kab_eq3=Pab_eq3*Cimu_obs_eq3'*inv(Cimu_obs_eq3*Pab_eq3*Cimu_obs_eq3'+Rimu_obs_eq3);
    Pab_eq3=Pab_eq3-Kab_eq3*Cimu_obs_eq3*Pab_eq3;
    dz=Kab_eq3*([imu_obs_eq2]-Cimu_obs_eq3*xest_ab_eq3);
    xest_ab_eq3=xest_ab_eq3+dz;
    
    if (flag)   %skip the first iteration because it is singular.
        flag=0; 
    else
        %Kalman equivalence
        Kab_eq=Pab_eq*Cst_eq'*inv(Cst_eq*Pab_eq*Cst_eq');
        Pab_eq=Pab_eq-Kab_eq*Cst_eq*Pab_eq;
        dz=Kab_eq*(-Cst_eq*xest_ab_eq);
        xest_ab_eq=xest_ab_eq+dz;
        
        %Kalman equivalence 1
        Kab_eq1=Pab_eq1*Cst_eq1'*inv(Cst_eq1*Pab_eq1*Cst_eq1');
        Pab_eq1=Pab_eq1-Kab_eq1*Cst_eq1*Pab_eq1;
        dz=Kab_eq1*(-Cst_eq1*xest_ab_eq1);
        xest_ab_eq1=xest_ab_eq1+dz;
    end
    
    debctr=debctr+1;
    debx(:,debctr)=xest;
    debxval(:,debctr)=[xn;x1;x2;x3;x4;x5];
    debp(:,debctr)=diag(P).^0.5;
    debxab(:,debctr)=xest_ab;
    debpab(:,debctr)=diag(Pab).^0.5;
    debxabcr(:,debctr)=xest_ab_cr;
    debpabcr(:,debctr)=diag(Pab_cr).^0.5;
    debxeq(:,debctr)=xest_ab_eq;
    debpeq(:,debctr)=diag(Pab_eq).^0.5;
    debxeq1(:,debctr)=xest_ab_eq1;
    debpeq1(:,debctr)=diag(Pab_eq1).^0.5;
    debxeq2(:,debctr)=xest_ab_eq2;
    debpeq2(:,debctr)=diag(Pab_eq2).^0.5;
    
    %Compute the states of next cycle
    xn=An*xn+N*u;
    x1=A1*x1+B1*randn(size(B1,1),1);
    x2=A2*x2+B2*randn(size(B2,1),1);
    x3=A3*x3+B3*randn(size(B3,1),1);
    x4=A4*x4+B4*randn(size(B4,1),1);
    x5=A5*x5+B5*randn(size(B5,1),1);
    
    %prediction for the next states
    xest=A*xest;
    xest(1:nst_n)=xest(1:nst_n)+N*imu_inp;
    P=A*P*A'+Q;
    
    xest_ab=Aab*xest_ab;
    xest_ab(1:nst_n)=xest_ab(1:nst_n)+N*imu_inpa;
    xest_ab((nst_n+1):(2*nst_n))=xest_ab((nst_n+1):(2*nst_n))+N*imu_inpb;
    Pab=Aab*Pab*Aab'+Qab;
    
    xest_ab_cr=Aab_cr*xest_ab_cr;
    xest_ab_cr(1:nst_n)=xest_ab_cr(1:nst_n)+N*imu_inpa;
    xest_ab_cr((nst_n+1):(2*nst_n))=xest_ab_cr((nst_n+1):(2*nst_n))+N*imu_inpb;
    xest_ab_cr=xest_ab_cr+N_cr*imu_obs;
    Pab_cr=Aab_cr*Pab_cr*Aab_cr'+Qab_cr;
    
    xest_ab_eq=Aab_eq*xest_ab_eq;
    xest_ab_eq(1:nst_n)=xest_ab_eq(1:nst_n)+N*imu_inpa;
    xest_ab_eq((nst_n+1):(2*nst_n))=xest_ab_eq((nst_n+1):(2*nst_n))+N*imu_inpb;
    Pab_eq=Aab_eq*Pab_eq*Aab_eq'+Qab_eq;
    
    xest_ab_eq1=Aab_eq1*xest_ab_eq1;
    xest_ab_eq1(1:nst_n)=xest_ab_eq1(1:nst_n)+N*imu_inpa;
    xest_ab_eq1((nst_n+1):(2*nst_n))=xest_ab_eq1((nst_n+1):(2*nst_n))+N*imu_inpb;
    Pab_eq1=Aab_eq1*Pab_eq1*Aab_eq1'+Qab_eq1;
    
    xest_ab_eq2=Aab_eq2*xest_ab_eq2;
    xest_ab_eq2(1:nst_n)=xest_ab_eq2(1:nst_n)+N*imu_inpa;
    xest_ab_eq2((nst_n+1):(2*nst_n))=xest_ab_eq2((nst_n+1):(2*nst_n))+N*imu_inpb;
    Pab_eq2=Aab_eq2*Pab_eq2*Aab_eq2'+Qab_eq2;
    
    xest_ab_eq3=Aab_eq3*xest_ab_eq3;
    xest_ab_eq3(1:nst_n)=xest_ab_eq3(1:nst_n)+N*imu_inpa;
    xest_ab_eq3((nst_n+1):(2*nst_n))=xest_ab_eq3((nst_n+1):(2*nst_n))+N*imu_inpb;
    Pab_eq3=Aab_eq3*Pab_eq3*Aab_eq3'+Qab_eq3;
    
    debxeq3(:,debctr)=xest_ab_eq3;
    debpeq3(:,debctr)=diag(Pab_eq3).^0.5;
end

figure(2);
yy=4; 
plot(debx(yy,1:debctr));
hold on;
plot(debxeq2(yy+nst_n,1:debctr));
plot(debxabcr(yy+nst_n,1:debctr));
plot(debxeq1(yy+nst_n,2:debctr),'r');
grid;
ylabel('Acc^2 bias Estimate (m/sec^2)');
xlabel('time (sec)');
legend('Sol 1','Sol 5a','Sol 5b', 'Sol 2');


figure(3);
yy=5; 
plot(debx(yy,1:debctr));
hold on;
plot(debxeq(yy+nst_n,1:debctr),'--r');
grid;
ylabel('Acc^3 bias Estimate (m/sec^2)');
xlabel('time (sec)');
legend('Sol 1','Sol 3');


figure(4);
yy=6; 
plot(debx(yy,1:debctr));
hold on;
plot(debxeq1(yy+nst_n,2:debctr),'k');
plot(debxeq3(yy+nst_n,1:debctr),'k');
grid;
ylabel('Acc^4 bias Estimate(m/sec^2)');
xlabel('time (sec)');
legend('Sol 1','Sol 2', 'Sol 4 (predicted)');
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
