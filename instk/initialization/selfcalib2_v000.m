%Sensor Refs=[use (0,1), val, obs noise, ini std]
%u is assumed to be deterministic unknown for each cycle
function [u, ximu, Pimu, Pu, Pux]=selfcalib2(imu_data, Aimu, Qimu, Cimu, Rimu, Pimu, M, sensor_ref)

DBG=1;

nu=size(M,2); %# of kinematic variables
nimu=size(Aimu,1);    %total # of imu states
ns=size(M,1); %# of sensors

%System model
A=Aimu;
Q=Qimu;
C=Cimu;
R=Rimu;
P=Pimu;

%External body acc/rot rate observations
yext=[];
for in=1:nu
    if (sensor_ref(in,1)==1)     %there is an external obs
        %Modify observation model
        M=[M;zeros(1,nu)];
        M(end,in)=1;
        C=[C;zeros(1,nimu)];
        R=diagmat_v000(sensor_ref(in,3)^2,R,[]);
        yext=[yext;sensor_ref(in,2)];
        ns=ns+1;
        if (sensor_ref(in,4)~=0)    %external obs have random bias
            A=diagmat_v000(1, A, []);
            Q=diagmat_v000(0, Q, []);
            P=diagmat_v000(sensor_ref(in,4)^2, P, []);
            C=[C zeros(ns,1)];
            C(end,end)=1;
            nimu=nimu+1;
        end
    end
end

%%Redundancy observation model
[T Mls]=cp_T_v000(M, R, 1); 
Cobs=T*C;
Robs=T*R*T';

if (DBG)
    debmls=zeros(nu,size(imu_data,2));
    debmls1=zeros(nu,size(imu_data,2));
    mx_a=M(1:size(imu_data,1),:);
    mx_b=inv(Rimu);
    Mls1=inv(mx_a'*mx_b*mx_a)*mx_a*mx_b;
    debop=zeros(nu,size(imu_data,2));
end

%%%%% Start the estimation routine
%Run the Kalman filter
nst=size(A,1);
xest=zeros(nst,1);
for in=1:size(imu_data,2)
    %%Filter
    y=T*[imu_data(:,in);yext];
    K=P*Cobs'*inv(Cobs*P*Cobs'+Robs);
    xest=xest+K*(y-Cobs*xest);
    P=(eye(nst)-K*Cobs)*P*(eye(nst)-K*Cobs)'+K*Robs*K';
    
    %%Predict
    xest=A*xest;
    P=A*P*A'+Q;
    
    if (DBG)
        debmls(:,in)=Mls*[imu_data(:,in);yext];
        debmls1(:,in)=Mls1*imu_data(:,in);
        debop(:,in)=Mls*([imu_data(:,in);yext]-C*xest);
    end
end

u=Mls*([imu_data(:,in);yext]-C*xest);
nimu=size(Aimu,1);
ximu=xest(1:nimu);
Pimu=P(1:nimu,1:nimu);
Pu=Mls*(C*P*C'+R)*Mls';
mx_a=Mls*C*P;
Pux=mx_a(:,1:nimu);

if (DBG)
    in=1;
    figure;
    plot(imu_data(in,:));
    hold on
    plot(debmls(in,:),'r');
    plot(debmls1(in,:),'k');
    plot(debop(in,:),'g');
end