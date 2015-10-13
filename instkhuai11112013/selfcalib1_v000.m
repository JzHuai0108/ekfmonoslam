%Sensor Refs=[use (0,1), val, std dev]
%Assumes u[k+1]=u[k];

function [u, ximu, Pimu, Pu, Pux]=selfcalib1(imu_data, Aimu, Qimu, Cimu, Rimu, Pimu, M, sensor_ref)

DBG=1;  %Debug

%under some conditions (guess which), all these operations are equivalent
%to the simple averaging... (u=mean(imu_data))

nu=size(M,2); %#of kinematic variables
nimu=size(Aimu,1);    %total # of states

if (DBG)
    ndat=size(imu_data,2);
    debmn=zeros(nu,ndat);
    debop=zeros(nimu+nu,ndat);
end

if (isempty(sensor_ref)) %not a compulsory argument
    sensor_ref=zeros(size(M,2),3);
end

%System model
A=Aimu;
Q=Qimu;
C=Cimu;
R=Rimu;
Pinv=inv(Pimu);
nst=nimu;   %number of observation

%%Seperate ref_det & adjust the overall model
ref_det=zeros(nu,1);
in_det=[];
in_sto=[];
for in=1:nu
    if (sensor_ref(in,1)==1 && sensor_ref(in,3)==0)     %perfectly known kinematic var.
        ref_det=ref_det+M(:,in)*sensor_ref(in,2);
        in_det=[in_det in];
    else   %kinematic var is not known exactly
        %Add it to the system model as a random constant with infinite
        %initial variance
        nst=nst+1;
        A=diagmat_v000(1,A,[]);
        Q=diagmat_v000(0,Q,[]);
        C=[C M(:,in)];
        Pinv=diagmat_v000(0,Pinv,[]);
        if (sensor_ref(in,1)==1 && sensor_ref(in,3)~=0)     %There is a noisy obs for kinematic var
            %Add it to the observation model
            C=[C;zeros(1,nst)];
            C(end,end)=1;
            R=diagmat_v000(sensor_ref(in,3)^2,R,[]);            
        end
        in_sto=[in_sto in];
    end
end

%%%%% Start the estimation routine
%%For the first data use the innovation form (Pinv is probably singular)
y=imu_data(:,1)-M*ref_det;
%xest=[dximu;zeros(nst-nimu),1)]
xest=zeros(nst,1);
Rinv=inv(R);
P=inv(Pinv+C'*Rinv*C);
xest=P*(Pinv*xest+C'*Rinv*y);

if (DBG)
    WLS=lsmat_v001(M,Rimu,0);
    mean_dat=WLS*y;
    debmn(:,1)=mean_dat;
    debop(:,1)=xest(1:nst);
end
    
%For the rest use standard KF (I just don't like the information form)
for in=2:size(imu_data,2) 
    %prediction
    xest=A*xest;
    P=A*P*A'+Q;
    
    %filter
    y=imu_data(:,in)-M*ref_det;
    K=P*C'*inv(C*P*C'+R);
    %K(1:nimu,:)=0;
    xest=xest+K*(y-C*xest);
    P=(eye(nst)-K*C)*P*(eye(nst)-K*C)'+K*R*K';
    
    if (DBG)
        mean_dat=((in-1)*mean_dat+(WLS*y))/in;
        debmn(:,in)=mean_dat;
        debop(:,in)=xest(1:nst);
    end
end

%Outputs (combine deterministic and stochastic results)
ximu=xest(1:nimu,1);
Pimu=P(1:nimu,1:nimu);
u=zeros(nu,1);
u(in_det,1)=sensor_ref(in_det,2);
u(in_sto,1)=xest(nimu+1:nst,:);
Pu=zeros(nu);
Pu(in_sto,in_sto)=P(nimu+1:nst,nimu+1:nst);
Pux=zeros(nu,nimu);
Pux(in_sto,:)=P(nimu+1:nst,1:nimu);


if (DBG)
    in=1;
    plot(imu_data(in,:))
    hold on;
    plot(debmn(in,:),'r')
    plot(debop(nimu+in,:),'g')
end


