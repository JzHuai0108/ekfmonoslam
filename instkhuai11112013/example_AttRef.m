%%NOTE: this example uses alingc_pl_v000 which assumes that pitch~=90.

clear;
fclose all;
dt=1/20;    %main frequency
motwin=0.7;     %sec. This must be smaller than the minimum time constant of imu error models.

motin=round(motwin/dt);
instat=round(0.5/dt);     %initial stationary period. 

%%Nominal nav values
gravity=9.808; %Gravity (assumed to be constant)
heading=10*pi/180;  %=magnetometer output. if there is no mag, use any arbitrary value.
heading_std=5/180*pi;   

%%interpolate the inputs to 20Hz
DIRNAME='C:\Benisil\PDANav\MotTest2\';
INPFILES=[[DIRNAME '12-17-42-02-acclog.bin ']; [DIRNAME '12-17-42-02-gyrolog.bin']];
tref=interp_imu_v001(INPFILES, [DIRNAME 'imu.bin'], [5 2 3 4 5;5 2 3 4 5], dt);

%%Read the input and calibrate
load('Sensors\calib_101008.mat');
imu=readbin_v000([DIRNAME 'imu.bin'],7);
Aaccinv=inv(Aacc);
Agyroinv=inv(Agyro);
CalA=[Aaccinv zeros(3,3);-Agyroinv*Bgyro*Aaccinv Agyroinv];
Calb=-[Aaccinv*bacc;Agyroinv*(Bgyro*Aaccinv*bacc+bgyro1+bgyro2)];
cimu=CalA*imu(2:7,:)+Calb*ones(1,size(imu,2));

%%Sensor Error model
load('Sensors\SenErrDef_20hz.mat');
[Aimu Bimu Cimu Dimu sPimu]=imu_modTI_v000(SenErrDef);
Qimu=Bimu*Bimu';
Rimu=Dimu*Dimu';

%%zero motion threshold (don't critize. I know it is not correct)
th=9*diag(Rimu)/2*4;

%Compute the initial attitude and self-calibrate
sensor_ref=[0 0 1e-6 1e-4;0 0 1e-6 1e-4;0 -gravity 1e-6 1e-1;1 0 1e-8 1e-5;1 0 1e-8 1e-5;1 0 1e-8 1e-5;];
[u, ximu, Pimu, Pu, Pux]=selfcalib1_v001(cimu(:,1:instat), Aimu, Qimu, Cimu, Rimu, inv(sPimu*sPimu'), eye(6), sensor_ref);
[Cpb E]=alingc_pl_v000(u(1:3), heading);
E=[E zeros(3,3)];
%initial attitude
Cbn=Cpb';

%Covariance matrix
nst=3+size(Aimu,1);
P=zeros(nst);
P(1:3,1:3)=E*Pu*E';
P(3,3)=P(3,3)+heading_std^2;
P(4:end,4:end)=Pimu;
P(1:3,4:end)=E*Pux;
P(4:end,1:3)=(E*Pux)';

%Accumulators
mctr=0;
macc=zeros(6,1);
pacc=zeros(6,1);

%Start navigation
debn=zeros(nst,size(cimu,2));
debP=zeros(nst,size(cimu,2));

for in=instat+1:size(cimu,2)
    %update the attitude
    sdat=cimu(:,in)-Cimu*ximu;
    
    rot=sdat(4:6)*dt;       %% Quess what can be inserted here?
    Cupd=rot2dcm_v000(rot);
    Cbn=Cbn*Cupd;

    %update the covariance
    Anav=eye(3);
    N=[zeros(3,3) -Cbn];
    STM=[Anav N*Cimu*dt; zeros(size(Aimu,1), 3) Aimu];
    Q=diagmat_v000(Qimu, N*Rimu*dt*N'*dt,1);
    P=STM*P*STM'+Q;
    
    %propagate the states
    ximu=Aimu*ximu;
    
    %zero motion detection
    %%TO DO: Implement a decent optimal motion detection algorithm!!!
    if (sum((cimu(:,in)-cimu(:,in-1)).^2<th,1)==6)
        macc=macc+sdat/motin;
        pacc=pacc+sdat.^2/motin;
        mctr=mctr+1;
    else
        mctr=0;
        macc=zeros(6,1);
        pacc=zeros(6,1);
    end
    
    %acc based orientation update
    if (mctr==motin)    %no instantenous motion detected in the last window
        if (1)  %check for the motion for the entire window //TBDL
            %No motion detected apply KF update
            heading=dcm2heading_v000(Cbn);
            [Cpb E]=alingc_pl_v000(macc(1:3), heading);
            
            %Observation
            mx_a=Cbn*Cpb;
            y=[mx_a(3,2);-mx_a(3,1);mx_a(2,1)];
            
            %Observation model
            H=zeros(3,nst);
            H(1,1)=-1;
            H(2,2)=-1;
            H(1:3,4:end)=-[E zeros(3,3)]*Cimu;
            
            R=Rimu(1:3,1:3)/motin;
            
            %Apply KF
            K=(P*H')/(H*P*H'+R);
            dz=K*(y);
            Cbn=(eye(3)+skew(dz(1:3)))*Cbn;
            ximu=ximu+dz(4:end);
            P=(eye(nst)-K*H)*P;
        end
        mctr=0;
        macc=zeros(6,1);
        pacc=zeros(6,1);
    end
    
    %%Appply magnometer updates here
    %
    %ipod does not have magnetometer and i do not have enough money to buy
    %iphone.
    %Waiting for android.
    
    debn(:,in)=[dcm2euler_v000(Cbn);ximu];
    debP(:,in)=diag(P).^0.5;
end

return;