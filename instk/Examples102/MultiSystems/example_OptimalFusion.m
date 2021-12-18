clear;

%Sensors
%%Acc set1
A1=diag([0.999987092709754 1]);
B1=diag([3.5596949329112e-006 1e-7]);
C1=[1 1];
sR1=(8.55478077790315e-007)^0.5;
sP1=diag([(4.9087e-007)^0.5 1e-3]);
[A1, Q1, R1]=dc2dc_v000(A1,B1*B1',sR1*sR1',1/100,2,1);
B1=Q1^0.5;
sR1=R1^0.5;


A2=A1;A3=A1;A4=A1;A5=A1;
B2=B1;B3=B1;B4=B1;B5=B1;
C2=C1;C3=C1;C4=C1;C5=C1;
sR2=sR1;sR3=sR1;sR4=sR1;sR5=sR1;
sP2=sP1;sP3=sP1;sP4=sP1;sP5=sP1;


%%Acc set2
A11=1;
B11=1e-7;%5.35244311625171e-008;
C11=1;
sR11=(8.55478077790315e-007)^0.5;
sP11=(6.8307e-007*20)^0.5;
[A11, Q11, R11]=dc2dc_v000(A11,B11*B11',sR11*sR11',1/100,2,1);
B11=Q11^0.5;
sR11=R11^0.5;

A12=A11;A13=A11;A14=A11;A15=A11;
B12=B11;B13=B11;B14=B11;B15=B11;
C12=C11;C13=C11;C14=C11;C15=C11;
sR12=sR11;sR13=sR11;sR14=sR11;sR15=sR11;
sP12=sP11;sP13=sP11;sP14=sP11;sP15=sP11;


%gyros
A6=0.999954464990211;
B6=7.50002159071526e-007*10;
C6=1;
sR6=(4.40958675010914e-008)^0.5;
vr_a=markov1st_v000(A6 ,B6,[], 2,100);
sP6=(vr_a(3))^0.5;
[A6, Q6, R6]=dc2dc_v000(A6,B6*B6',sR6*sR6',1/100,2,1);
B6=Q6^0.5;
sR6=R6^0.5;

A7=A6;B7=B6;C7=C6;sR7=sR6;sP7=sP6;

% A7=0.999977286299473;
% B7=2.06907835784507e-006;
% C7=1;
% sR7=(3.74167215805145e-007)^0.5;
% sP7=(9.4241e-008)^0.5;
% [A7, Q7, R7]=dc2dc_v000(A7,B7*B7',sR7*sR7',1/100,2,1);
% B7=Q7^0.5;
% sR7=R7^0.5;

MA1=[cos(0) sin(0) 0;
    cos(0) sin(0) 0;
    cos(-pi/5) sin(-pi/5) 0;
    cos(-pi/2) sin(-pi/2) 0;
    cos(-pi/2) sin(-pi/2) 0];


MG1=[0 0 1;0 0 1];

MA2=[cos(pi/10) sin(pi/10) 0;
    cos(pi/8) sin(pi/8) 0;
    cos(pi/6) sin(pi/6) 0;
    cos(pi/4) sin(pi/4) 0;
    cos(pi/2) sin(pi/2) 0];

%%% Main model
%M=[MA1;MG];
M=[MA1;MG1;MA2];

%imu observation structure
Rimu=sR1*sR1';
Rimu=diagmat_v000(sR2*sR2',Rimu,1);
Rimu=diagmat_v000(sR3*sR3',Rimu,1);
Rimu=diagmat_v000(sR4*sR4',Rimu,1);
Rimu=diagmat_v000(sR5*sR5',Rimu,1);
Rimu=diagmat_v000(sR6*sR6',Rimu,1);
Rimu=diagmat_v000(sR7*sR7',Rimu,1);
Rimu=diagmat_v000(sR11*sR11',Rimu,1);
Rimu=diagmat_v000(sR12*sR12',Rimu,1);
Rimu=diagmat_v000(sR13*sR13',Rimu,1);
Rimu=diagmat_v000(sR14*sR14',Rimu,1);
Rimu=diagmat_v000(sR15*sR15',Rimu,1);

Cimu=C1;
Cimu=diagmat_v000(C2,Cimu,1);
Cimu=diagmat_v000(C3,Cimu,1);
Cimu=diagmat_v000(C4,Cimu,1);
Cimu=diagmat_v000(C5,Cimu,1);
Cimu=diagmat_v000(C6,Cimu,1);
Cimu=diagmat_v000(C7,Cimu,1);
Cimu=diagmat_v000(C11,Cimu,1);
Cimu=diagmat_v000(C12,Cimu,1);
Cimu=diagmat_v000(C13,Cimu,1);
Cimu=diagmat_v000(C14,Cimu,1);
Cimu=diagmat_v000(C15,Cimu,1);

[T Mls]=cp_Tparam_v000(M,Rimu);
Cimu_obs=[T*Cimu];
Rimu_obs=T*Rimu*T';
Rimu_inp=Mls*Rimu*Mls';


A=A1;
A=diagmat_v000(A2,A,1);
A=diagmat_v000(A3,A,1);
A=diagmat_v000(A4,A,1);
A=diagmat_v000(A5,A,1);
A=diagmat_v000(A6,A,1);
A=diagmat_v000(A7,A,1);
A=diagmat_v000(A11,A,1);
A=diagmat_v000(A12,A,1);
A=diagmat_v000(A13,A,1);
A=diagmat_v000(A14,A,1);
A=diagmat_v000(A15,A,1);

Q=B1*B1';
Q=diagmat_v000(B2*B2',Q,1);
Q=diagmat_v000(B3*B3',Q,1);
Q=diagmat_v000(B4*B4',Q,1);
Q=diagmat_v000(B5*B5',Q,1);
Q=diagmat_v000(B6*B6',Q,1);
Q=diagmat_v000(B7*B7',Q,1);
Q=diagmat_v000(B11*B11',Q,1);
Q=diagmat_v000(B12*B12',Q,1);
Q=diagmat_v000(B13*B13',Q,1);
Q=diagmat_v000(B14*B14',Q,1);
Q=diagmat_v000(B15*B15',Q,1);


P=sP1*sP1';
P=diagmat_v000(sP2*sP2',P,1);
P=diagmat_v000(sP3*sP3',P,1);
P=diagmat_v000(sP4*sP4',P,1);
P=diagmat_v000(sP5*sP5',P,1);
P=diagmat_v000(sP6*sP6',P,1);
P=diagmat_v000(sP7*sP7',P,1);
P=diagmat_v000(sP11*sP11',P,1);
P=diagmat_v000(sP12*sP12',P,1);
P=diagmat_v000(sP13*sP13',P,1);
P=diagmat_v000(sP14*sP14',P,1);
P=diagmat_v000(sP15*sP15',P,1);

Pwls=P;

% %%%model for the equivalent systems
%Equivalent models for identical systems
[Mls_e1, Ae1, Be1, Ce1, sRe1, sPe1]=cp_redsys_v000(MA1(:,1:2), A1, B1, C1, sR1, sP1);
[Mls_e2, Ae2, Be2, Ce2, sRe2, sPe2]=cp_redsys_v000(MA2(:,1:2), A11, B11, C11, sR11, sP11);
[Mls_e3, Ae3, Be3, Ce3, sRe3, sPe3]=cp_redsys_v000(MG1(:,3), A6, B6, C6, sR6, sP6);

%imu observation structure
Rimue=sRe1*sRe1';
Rimue=diagmat_v000(sR6*sR6',Rimue,1);
Rimue=diagmat_v000(sR7*sR7',Rimue,1);
Rimue=diagmat_v000(sRe2*sRe2',Rimue,1);

Cimue=Ce1;
Cimue=diagmat_v000(C6,Cimue,1);
Cimue=diagmat_v000(C7,Cimue,1);
Cimue=diagmat_v000(Ce2,Cimue,1);

Me=[1 0 0;0 1 0;MG1;1 0 0;0 1 0];
[Te Mlse]=cp_Tparam_v000(Me,Rimue);
Cimu_obse=[Te*Cimue];
Rimu_obse=Te*Rimue*Te';
Rimu_inpe=Mlse*Rimue*Mlse';

%single system model
Ae=Ae1;
Ae=diagmat_v000(A6,Ae,1);
Ae=diagmat_v000(A7,Ae,1);
Ae=diagmat_v000(Ae2,Ae,1);

Qe=Be1*Be1';
Qe=diagmat_v000(B6*B6',Qe,1);
Qe=diagmat_v000(B7*B7',Qe,1);
Qe=diagmat_v000(Be2*Be2',Qe,1);

Pe=sPe1*sPe1';
Pe=diagmat_v000(sP6*sP6',Pe,1);
Pe=diagmat_v000(sP7*sP7',Pe,1);
Pe=diagmat_v000(sPe2*sPe2',Pe,1);


%%%%%%%%%%%%%%%%%%
nit=60*100;
nst=size(A,1);
nste=size(Ae,1);

x=(P^0.5)*randn(nst,1);
xest=zeros(size(A,1),1);
xeste=zeros(size(Ae,1),1);

sRimu=Rimu^0.5;
nsen=size(sRimu,1);
sQ=Q^0.5;

%debug variables
debctr=0;
debx=zeros(nst,nit);
debp=zeros(nst,nit);
debxe=zeros(nste,nit);
debpe=zeros(nste,nit);
debuopt=zeros(3,nit);
debuopte=zeros(3,nit);
debuwls=zeros(3,nit);
debpopt=zeros(3,nit);
debpopte=zeros(3,nit);
debpwls=zeros(3,nit);
debeq=zeros(5,nit);


% debctr=debctr+1;
% debx(:,debctr)=xest;
% debp(:,debctr)=diag(P).^0.5;
% debxe(:,debctr)=xeste;
% debpe(:,debctr)=diag(Pe).^0.5;

%Start main loop
for in=1:nit
    %observations
    u=[0;0;0];
    
    y=M*u+Cimu*x+sRimu*randn(nsen,1);
    imu_inp=Mls*y;
    imu_obs=T*y;
    
    ye=[Mls_e1*y(1:5);y(6:7);Mls_e2*y(8:end)];
    imu_inpe=Mlse*ye;
    imu_obse=Te*ye;
    
    %Kalman
    K=P*Cimu_obs'*inv(Cimu_obs*P*Cimu_obs'+Rimu_obs);
    P=P-K*Cimu_obs*P;
    dz=K*(imu_obs-Cimu_obs*xest);
    xest=xest+dz;
    
    Ke=Pe*Cimu_obse'*inv(Cimu_obse*Pe*Cimu_obse'+Rimu_obse);
    Pe=Pe-Ke*Cimu_obse*Pe;
    dze=Ke*(imu_obse-Cimu_obse*xeste);
    xeste=xeste+dze;
    
    mx_a=Mls*Cimu*P*Cimu'*Mls';
    mx_b=Mlse*Cimue*Pe*Cimue'*Mlse';
    mx_c=Mls*Cimu*Pwls*Cimu'*Mls';
    
    debctr=debctr+1;
    debx(:,debctr)=xest;
    debp(:,debctr)=diag(P).^0.5;
    debxe(:,debctr)=xeste;
    debpe(:,debctr)=diag(Pe).^0.5;
    debuopt(:,debctr)=Mls*y-Mls*Cimu*xest;
    debuopte(:,debctr)=Mlse*ye-Mlse*Cimue*xeste;
    debuwls(:,debctr)=Mls*(y);
    debpopt(:,debctr)=diag(mx_a+Rimu_inp).^0.5;
    debpopte(:,debctr)=diag(mx_b+Rimu_inpe).^0.5;
    debpwls(:,debctr)=diag(mx_c+Rimu_inp).^0.5;
    debeq(1:2,debctr)=Mls_e1*y(1:5);
    debeq(3,debctr)=Mls_e3*y(6:7);
    debeq(4:5,debctr)=Mls_e2*y(8:end);
    
    %Compute the states of next cycle
    x=A*x+sQ*randn(nst,1);
           
    %prediction for the next states
    xest=A*xest;
    P=A*P*A'+Q;
    Pwls=A*Pwls*A'+Q;
    
    xeste=Ae*xeste;
    Pe=Ae*Pe*Ae'+Qe;
    
end

figure(10);
yy=1;
tt=[1:10:debctr]/60;
plot(tt,debuopte(yy,1:10:debctr));
hold on;
plot(tt,debuopt(yy,1:10:debctr));
plot(tt,debuwls(yy,1:10:debctr),'k');
plot(tt,debeq(0+yy,1:10:debctr),'r');
plot(tt,debeq(3+yy,1:10:debctr),'m');

plot(tt,debpopte(yy,1:10:debctr),'--');
plot(tt,debpopt(yy,1:10:debctr),'--');
plot(tt,debpwls(yy,1:10:debctr),'--k');
grid;
xlabel('Time (min)');
ylabel('X-Acceleration Error(m/sec^2)');
%legend('Opt. S. (u^{opt})', 'Red. Opt. S. (u^{red})','WLS S. (u^{wls})','Set1 Opt. (u^{A1})','Set2 Opt. (u^{A2})','Opt. SD.', 'Red. Opt. SD.','WLS SD');
legend('Opt. S. (u^{opt})', 'Red. Opt. S. (u^{red})','WLS S. (u^{wls})','Set1 Opt. (u^{A1})','Set2 Opt. (u^{A2})');

figure(11);
yy=2;
tt=[1:10:debctr]/60;
plot(tt,debuopte(yy,1:10:debctr));
hold on;
plot(tt,debuopt(yy,1:10:debctr));
plot(tt,debuwls(yy,1:10:debctr),'k');
plot(tt,debeq(0+yy,1:10:debctr),'r');
plot(tt,debeq(3+yy,1:10:debctr),'m');

plot(tt,debpopte(yy,1:10:debctr),'--');
plot(tt,debpopt(yy,1:10:debctr),'--');
plot(tt,debpwls(yy,1:10:debctr),'--k');
grid;
xlabel('Time (min)');
ylabel('Y-Acceleration Error(m/sec^2)');
%legend('Opt. S. (u^{opt})', 'Red. Opt. S. (u^{red})','WLS S. (u^{wls})','Set1 Opt. (u^{A1})','Set2 Opt. (u^{A2})','Opt. SD.', 'Red. Opt. SD.','WLS SD');
legend('Opt. S. (u^{opt})', 'Red. Opt. S. (u^{red})','WLS S. (u^{wls})','Set1 Opt. (u^{A1})','Set2 Opt. (u^{A2})');

return;
