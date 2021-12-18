%%Forward solution of RTS.
%In this solution I directly save the
%predicted and updated states so that I can apply RTS formula as it is.
%This is the simplest smoother implementation one can choose. In the
%forward path I only add 2 lines to record the updated and predicted values
%(see the lines starting with "smoother stuff"). That is the only change
%required in the forward path of this type of RTS implementation. The rest is handled in the
%backwards part.
%PS:You do not have to save the any "updated solution" when there is no
%update. In fact, it would be much better to record predicted and updated solution
%seperetaly. In the following code, I record updated solution in each cycle
%(regardless of whether or not there is a real update) just to make it
%easier to follow and understand the code by eye. In a real implementation
%of course you must record minimum number of variables at each cycle.

fclose all;
clear;

%Time configuration
TIME.interval=[20 450];    %start/stop times
TIME.dt=1/400;  %sampling time (required only to seek faster. dt is actually computed from the file)
TIME.salign=5;   %static alignment duration for the beginning
TIME.zuptInterval=[21 30];    %Time intervals to apply zupt
TIME.zuptDt=0.2;  %minimum time between 2 zupt updates
TIME.maxCovStep=0.1; %maximum covariance propagation step

%Observation params
OBS.odoR=0.2^2;  %(m/s)^2
OBS.nonHoloR=0.1^2; %(m/s)^2
OBS.gpsR=0.3^2'; %(m)^2
OBS.zuptR=0.1^2;
OBS.LA=[0.2;0;0];

IMUTYPE=1;

INDIR='/InputDirectory/I/want/to/live/in/Shanghai/';
OUTDIR='/OutputDirectory/I/want/to/work/in/Sydney/'; %Yet, life force me to dwell in Ankara.

%open files
fimu=fopen([INDIR 'sensordata.bin'],'r');
imudata=fread(fimu,7,'double');
sr_a=round((TIME.interval(1)-imudata(1))/TIME.dt);
if sr_a>0
    fseek(fimu, sr_a*7*8, 0);
end

fodo=fopen([INDIR 'odo.txt'],'rt');
fgps=fopen([INDIR 'gps.txt'],'rt');

ffilres=fopen([OUTDIR 'navresult.bin'],'w');
fsmobuf=fopen([OUTDIR 'smobuffer.bin'],'w');

%Initial PVA
iniLlh=[39.88864*pi/180;32.78002*pi/180;1112];
Llh=iniLlh;
Vn=[0;0;0];
qbn=[1;0;0;0]; %Unknown:To be determined during coarse alignment

%Load IMU error definitions
IMU_ERRDEF=imu_err_defs_v000(IMUTYPE); %needed just to initialize P. The rest is handled within sys model

%Initial Covariance matrix
nst=21;
P=zeros(nst);
[Rn, Re, ~, sL, cL, WIE_E]=geoparam_v000(Llh);
Rappx=sqrt((Rn+Llh(3))*(Re+Llh(3)));
P(1:2,1:2)=eye(2)*(5/Rappx)^2; %position errors in meter
P(3,3)=5^2; %position errors in meter
P(4:6,4:6)=eye(3)*0.05^2; %vel error in m/s
P(7:9,7:9)=eye(3); %attitude errors in radian (will be determined later during coarse alignment)
P(10:12,10:12)=eye(3)*IMU_ERRDEF.acc_bias_var;
P(13:15,13:15)=eye(3)*IMU_ERRDEF.gyro_bias_var;
P(16:18,16:18)=eye(3)*IMU_ERRDEF.acc_scale_var;
P(19:21,19:21)=eye(3)*IMU_ERRDEF.gyro_scale_var;


%%%%%Static Alignment
%%Step I:Coarse align
coarse_dur=TIME.salign;
%Compute the mean of imu data during the coarse alignment period
datbuf=fread(fimu, [7 round(coarse_dur/TIME.dt)], 'double');
meandat=mean(datbuf(2:end,:),2)/TIME.dt;
preimutime=datbuf(1,end);   %the last imu time
clear datbuf;

%Compute the alignment between the geodetic and the body frame
[~, ~, wie_n g]=geoparam_v001(2, [cos(Llh(1));0;-sin(Llh(1))], Llh(3), zeros(3,1));
[Cnb E err]=align_opt_v000([meandat(4:6), meandat(1:3)],[[0;0;-g], wie_n] ,[1 0]); %Note:E defined for errors on imu. Therefore, the order of vectors is important
qbn=dcm2quat_v000(Cnb');

%%Step II:Determine the attitude covariance
%%Note I will ignore the initial cross correlation between the imu states
%%and the initial attitude. I will use convervative correlation values and
%%expect the initial ZUPT to take care of cross correlations.
P(7,7)=(IMU_ERRDEF.acc_vrw/10+IMU_ERRDEF.acc_bias_var)/9.81^2; %assume only a 1 sec average. 
P(8,8)=(IMU_ERRDEF.acc_vrw/10+IMU_ERRDEF.acc_bias_var)/9.81^2;
P(9,9)=(IMU_ERRDEF.gyro_arw/10+IMU_ERRDEF.gyro_bias_var)/(7.2e-5*cos(Llh(1)))^2;
P(9,9)=P(9,9)+(3*pi/180)^2; %forget about the flicker and use a random value for the initial bias error.


%%Step III:Self calibrate based on zero-motion
%TODO Later


%%%%%Discard all observation data before the current time
ododata=fscanf(fodo,'%f%f',[1,2]);
while (ododata(1)<=preimutime)
    ododata=fscanf(fodo,'%f%f',[1,2]);
end
gpsdata=fscanf(fgps,'%f%f',[1,7]);
while (gpsdata(1)<=preimutime)
    gpsdata=fscanf(fgps,'%f%f',[1,7]);
end


%%%%%Start the main INS
%prepare for the loop
covupt_time=preimutime;  %the time that we last updated the covariance
zupt_time=preimutime;   %the time that we last applied a zupt
imuaccum=zeros(6,1);
imuerrors=zeros(12,1);
gpsctr=0;

%I will accumulate 4 imu samples before processing.
imudata=fread(fimu, [7,4], 'double'); %imudata must be in increment form (not rate). If the file contains rate data, correct it here.
curimutime=imudata(1,end);

% ododata(1)=inf;
% TIME.zuptInterval=[inf 0;inf 0];

while (~feof(fimu))
    %Write result to the files   
    vr_a=posdiff_v000(Llh(1:2)', Llh(3), iniLlh);
    vr_c=quat2dcm_v000(qbn);
    fwrite(ffilres,[preimutime;vr_a';Vn;dcm2euler_v000(vr_c)*180/pi;imuerrors;diag(P)],'double');
    
    %downsample the data by accumulating it. (here I will downsample by a
    %factor of 4)
    angleinc=sum(imudata(2:4,:),2);
    velinc=sum(imudata(5:7,:),2);
    dt=curimutime-preimutime;
    
    %Calibrate the imu
    angleinc=angleinc-(imuerrors(4:6)*dt+diag(angleinc)*imuerrors(10:12)/1000);
    velinc=velinc-(imuerrors(1:3)*dt+diag(velinc)*imuerrors(7:9)/1000);
    
    %imu data accumulator for the covariance update
    imuaccum=imuaccum+[angleinc;velinc];
    
    %run strapdown with angle and velocity increments
    [qbn, Vn, Llh]=strapdown_ned_quat_v000(qbn, Vn, Llh, velinc, angleinc, dt);
    
    %Check if we are going to apply ZUPT.
    zupt_upd=0;
    if (curimutime-zupt_time>=TIME.zuptDt)
        if (curimutime>TIME.zuptInterval(1,1) && curimutime<TIME.zuptInterval(1,2)) %Note:If you have more than 1 zupt interval, use a for loop here to iterate all over
            zupt_upd=1;
        end
    end
    
    
    %Update the covariance
    smo_write=0;
    if ((curimutime-covupt_time)>=TIME.maxCovStep || curimutime>ododata(1) || curimutime>gpsdata(1) || zupt_upd)
        %propagate the covariance
        covdt=curimutime-covupt_time;
        [STM, Qd]=sys_ned_dcm_v001(Llh, Vn, qbn, imuaccum(4:6)/covdt, imuaccum(1:3)/covdt, covdt, IMUTYPE);
        P=STM*P*STM'+Qd;
        
        %Propagate the imu error states
        imuerrors=STM(10:end,10:end)*imuerrors;
        
        covupt_time=curimutime;
        imuaccum=zeros(6,1);
        
        %%Smoother Stuff
        %Record the predicted solutions
        predicted=[curimutime;Llh;Vn;qbn;imuerrors;mat2vec_v000(P);STM(:)];
        smo_write=1;
    end
    
    
    %Apply zupt
    if (zupt_upd)
        inno=Vn;
        H=[zeros(3) eye(3) zeros(3,15)];
        R=eye(3)*OBS.zuptR;

        %Kalman
        klmn_scrpt;

        zupt_time=curimutime;
    end


    %%Odometer updates
    if (curimutime>ododata(1))  %we have an odometer data
        Cnb=quat2dcm_v000(qbn)';
        inno=[ododata(2);0;0]-Cnb*Vn;
        H=[zeros(3) -Cnb Cnb*skew(Vn) zeros(3,12)];
        R=diag([OBS.odoR, OBS.nonHoloR, OBS.nonHoloR]);

        %Kalman
        klmn_scrpt;

        %Read the next odometer data
        ododata=fscanf(fodo,'%f%f',[1,2]);
        if (isempty(ododata))
            ododata(1)=inf;
        end
    end


    %%GPS updates
    if (curimutime>gpsdata(1))
        vr_a=posdiff_v000(Llh(1:2)', Llh(3), [gpsdata(2:3)*pi/180 gpsdata(4)]);
        disp(vr_a');

        lever=diag([1/Rappx, 1/Rappx,1])*quatrot_v000(qbn,OBS.LA,0); %Thx Jianzhu Huai
        inno=Llh-[gpsdata(2:3)*pi/180 gpsdata(4)]'-lever;
        H=[eye(3) zeros(3,18)]; %Completely ignore the effect of Cbn errors on level-arm compensation. If you want, use H(:,7:9) = skew(OBS.LA);
        R=diag([1/Rappx^2, 1/Rappx^2,1])*OBS.gpsR;

        %Kalman
        klmn_scrpt;

        %Read the new gps data
        gpsdata=fscanf(fgps,'%f%f',[1,7]);
        if (isempty(gpsdata))
            gpsdata(1)=inf;
        end

        gpsctr=gpsctr+1;
    end  
    
    %%Smoother stuff
    if (smo_write)
        updated=[curimutime;Llh;Vn;qbn;imuerrors;mat2vec_v000(P)]; %Note that if there is no update then updated==predicted. I simply ignore this fact that to reduce the number of code lines.
        fwrite(fsmobuf,[predicted;updated],'double');
        if (gpsctr==2)  %DELETE ME
            break;
        end
    end

    %Read the next imu data
    preimutime=curimutime;
    imudata=fread(fimu, [7,4], 'double');
    if (isempty(imudata))
        break;
    end
    curimutime=imudata(1,end);
end

%Last filtered solution
vr_a=posdiff_v000(Llh(1:2)', Llh(3), iniLlh);
vr_c=quat2dcm_v000(qbn);
fwrite(ffilres,[curimutime;vr_a';Vn;dcm2euler_v000(vr_c)*180/pi;imuerrors;diag(P)],'double');

fclose all;
