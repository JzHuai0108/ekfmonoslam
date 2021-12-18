fclose all;
clear;

%I/O Files
DIRPATH='\Tight\sequential\data\';
FIMU=fopen([DIRPATH 'IMU.txt'],'r'); %%Assumed to be sampled at stricly 400Hz
FGPS=fopen([DIRPATH 'GPS.txt'],'r'); %%In order not to deal with program control stutucture it is assumed that GPS data is available at strictly 1Hz.
FCOV=fopen([DIRPATH 'Cov.txt'],'wt');
FOUT=fopen([DIRPATH 'Output.txt'],'wt');

%Load System parameters
[SenErrDef, ClkErrDef, IniErrDef, ObsErrDef, ini_pva, dt]=sys_def_v000();

%initial values
pos_e=ecef2geo_v000(ini_pva(:,1),1);
Cne=pos2Cne_v000(ini_pva(1,1),ini_pva(2,1));
Cbn=euler2dcm_v000(ini_pva(:,3));
Cbe=Cne*Cbn;
vel_e=Cne*ini_pva(:,2);

%Sensor error models
[Aimu, Bimu, Cimu, Dimu, sP0imu]=imu_modTI_v000(SenErrDef);
Qimu=Bimu*Bimu';
Rimu=Dimu*Dimu';
%Clock Error model
[Aclk, Bclk, Cclk, Dclk, sP0clk]=imu_modTI_v000(ClkErrDef);

%%In fact this part should be performed inside the propagate() routine.
%However, as I fixed the the propagation period, I simply discretize it with
%the known dt;
mx_a = 40*dt*[-Aclk,Bclk*Bclk';zeros(size(Aclk)),Aclk'];   
mx_b = expm(mx_a);
STMclk = mx_b(3:4,3:4)';
Qclk = STMclk*mx_b(1:2,3:4);



%Initial covariance
nst=17;
P=zeros(nst);
P(1:3,1:3)=IniErrDef.pos_sP*IniErrDef.pos_sP';
P(4:6,4:6)=IniErrDef.vel_sP*IniErrDef.vel_sP';
P(7:9,7:9)=IniErrDef.att_sP*IniErrDef.att_sP';
P(10:15,10:15)=sP0imu*sP0imu';  %IMU errors
P(16:17,16:17)=sP0clk*sP0clk';  %Clock Errors

%Augmented states
Ximu=zeros(size(Aimu,1),1);
Xclk=zeros(size(Aclk,1),1);

pr_ctr=0;
imu_inc=zeros(6,1);
%Read the first IMU data
imu_data=fscanf(FIMU, '%f', 6);	%order=[Gyro;Acc]
while(~feof(FIMU))
    acc=imu_data(4:6)-Ximu(1:3);    %Correct them with error estimates
    gyro=imu_data(1:3)-Ximu(4:6);
    
    %Total IMU increments (used for low rate propagation)
    imu_inc=imu_inc+[acc;gyro]*dt;
    
    %INS calculations
    [Cbe, vel_e, pos_e]=strapdown_ecef_dcm_v000(Cbe, vel_e, pos_e, acc, gyro, dt);
    pr_ctr=pr_ctr+1;
    
    %Propagate the covariance and states at 10Hz
    if (~mod(pr_ctr,40))
        %Navigation model
        [STMNav QNav]=mdl_ecef_dcm_v000(Cbe, imu_inc(1:3)/0.1, Aimu, Qimu, Cimu, Rimu, dt*40);
        
        %Augment the clock model to nav model
        STM=diagmat_v000(STMclk,STMNav,1);
        Q=diagmat_v000(Qclk,QNav,1);
        
        %Propagate Covariance
        P=STM*P*STM'+Q;
        
        %Propagate IMU and Clock States
        Ximu=STMNav(10:15,10:15)*Ximu;
        Xclk=STMclk*Xclk;
        
        imu_inc=zeros(6,1);
    end
    
    %%Record the outputs
    if (~mod(pr_ctr,40))
        Llh=ecef2geo_v000(pos_e,0);
        Cne=pos2Cne_v000(Llh(1),Llh(2));
        vel_n=Cne'*vel_e;
        Cbn=Cne'*Cbe;
        eul=dcm2euler_v000(Cbn);
        fprintf(FOUT, '%12.12f ', [pr_ctr;Llh;vel_n;eul;Ximu;Xclk]);
        fprintf(FOUT, '\n');
        fprintf(FCOV, '%12.12f', [diag(P).^0.5]);
        fprintf(FCOV, '\n');
    end
        
    
    
    %Apply the GPS Pseudorange at 1Hz
    if (~mod(pr_ctr,400))
        gps_data=get_gps_data(FGPS);    %gps_data(i,:)=[sat_prn,pseudo distance,sat_x,sat_y,sat_z];
        
        %Start Sequential Processing of PseudoRanges
        dx=zeros(nst,1);
        for (in=1:size(gps_data,1))
            %Compute the innovation
            nom_range=norm(gps_data(in,3:5)'-pos_e);
            inno=gps_data(in,2)-(nom_range+Cclk*Xclk);
            
            %Observation Matrix
            H=zeros(1,nst);
            H(1,1:3)=(gps_data(in,3:5)-pos_e')/nom_range;
            H(1,16:17)=Cclk;
            
            Robs=ObsErrDef.sR*ObsErrDef.sR';
            
            %Update
            K=(P*H')/(H*P*H'+Robs);
            P=(eye(nst)-K*H)*P;
            dx=dx+K*(inno-H*dx);
        end
        
        %Correct states
        pos_e=pos_e-dx(1:3);
        vel_e=vel_e-dx(4:6);
        mx_a=eye(3)+skew(dx(7:9));
        Cbe=mx_a*Cbe;
        Ximu=Ximu+dx(10:15);
        Xclk=Xclk+dx(16:17);
    end
    
    %Read the next imu data from file
    imu_data=fscanf(FIMU, '%f', 6);
end

fclose(FIMU);
fclose(FGPS);
fclose(FCOV);
fclose(FOUT);

return;

%%plot the results
load([DIRPATH  'output.txt']);
plot(output(:,2),output(:,3));

        
            
        