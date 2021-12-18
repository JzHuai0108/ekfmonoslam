%Sensor Refs=[use (0,1), val, obs noise std, ini std]
%u is assumed to be a random constant (with infinite initial variance)
function [u, ximu, Pimu, Pu, Pux]=selfcalib1(imu_data, Aimu, Qimu, Cimu, Rimu, PimuInv, M, sensor_ref)

DBG=0;  %Debug

nu=size(M,2); %# of kinematic variables
ns=size(M,1); %# of sensors
nimu=size(Aimu,1);    %total # of imu states

if (isempty(sensor_ref)) %not a compulsory argument
    sensor_ref=zeros(size(M,2),4);
end

%System model (augment with kv)
A=diagmat_v000(Aimu, eye(nu),[]);
Q=diagmat_v000(Qimu, zeros(nu), []);
C=[M Cimu];
R=Rimu;
Pinv=diagmat_v000(PimuInv, zeros(nu), []);

%External body acc/rot rate observations
next=0;
Cext=[];
Rext=[];
yext=[];
for in=1:nu
    if (sensor_ref(in,1)==1)     %there is an external obs
        next=next+1;
        Cext=[Cext;zeros(1,size(A,1))];
        Cext(end,in)=1;
        Rext=diagmat_v000(sensor_ref(in,3)^2,Rext,[]);
        yext=[yext;sensor_ref(in,2)];
        if (sensor_ref(in,4)~=0)    %external obs have random bias
            A=diagmat_v000(1, A, []);
            Q=diagmat_v000(0, Q, []);
            Pinv=diagmat_v000(1/sensor_ref(in,4)^2, Pinv, []);
            C=[C zeros(ns,1)];
            Cext=[Cext zeros(next,1)];
            Cext(end,end)=1;
        end
    end
end
nst=size(A,1);

%%%%% Start the estimation routine
%%Sensor outputs
%%For the first sensor data use the innovation form (Pinv is probably singular)
y=imu_data(:,1);
Rinv=inv(R);
P=inv(Pinv+C'*Rinv*C);
xest=P*(C'*Rinv*y);
%%External obs
if (next>0)
    K=P*Cext'*inv(Cext*P*Cext'+Rext);
    xest=xest+K*(yext-Cext*xest);
    P=(eye(nst)-K*Cext)*P*(eye(nst)-K*Cext)'+K*Rext*K';
end


if (DBG)
    ndat=size(imu_data,2);
    debmn=zeros(nu,ndat);
    debop=zeros(size(A,1),ndat);
    
    WLS=lsmat_v001(M,R,0);
    mean_dat=WLS*y;
    debmn(:,1)=mean_dat;
    debop(:,1)=xest;
end

%For the rest use standard KF
for in=2:size(imu_data,2) 
    %%%prediction
    xest=A*xest;
    P=A*P*A'+Q;
    
    %%%filter
    %%Sensor data
    y=imu_data(:,in);
    K=P*C'*inv(C*P*C'+R);
    %K(1:nimu,:)=0;
    xest=xest+K*(y-C*xest);
    P=(eye(nst)-K*C)*P*(eye(nst)-K*C)'+K*R*K';
    %%External obs
    if (next>0)
        K=P*Cext'*inv(Cext*P*Cext'+Rext);
        xest=xest+K*(yext-Cext*xest);
        P=(eye(nst)-K*Cext)*P*(eye(nst)-K*Cext)'+K*Rext*K';
    end
    
    if (DBG)
        mean_dat=((in-1)*mean_dat+(WLS*y))/in;
        debmn(:,in)=mean_dat;
        debop(:,in)=xest;
    end
end

%Outputs (combine deterministic and stochastic results)
ximu=xest(nu+1:nu+nimu,1);
Pimu=P(nu+1:nu+nimu,nu+1:nu+nimu);
u=xest(1:nu,:);
Pu=P(1:nu,1:nu);
Pux=P(1:nu,nu+1:nu+nimu);


if (DBG)
    in=1;
    figure;
    plot(imu_data(in,:))
    hold on;
    plot(debmn(in,:),'r')
    plot(debop(in,:),'g')
end


