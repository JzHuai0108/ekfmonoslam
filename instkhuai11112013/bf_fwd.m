fclose all;
clear;

%Time configuration
TIME.interval=[10 2500];    %start/stop
TIME.dt=1/400;  %sampling time (required only to seek faster. dt is actually computed from the file)
TIME.salign=25;   %static alignment duration for the beginning
TIME.zuptInterval=[40 70];    %Time intervals to apply zupt
TIME.zuptDt=0.5;  %minimum time between 2 zupt updates
TIME.maxCovStep=0.1; %maximum covariance propagation step

%Observation params
OBS.odoR=0.2^2;  %(m/s)^2
OBS.useNonHolo=1;   %0:Don't use, 1:use
OBS.nonHoloR=0.1^2; %(m/s)^2
OBS.gpsR=0.1^2'; %(m)^2
OBS.zuptR=0.1^2;
OBS.LA=[0.2;0;0];

INDIR='/home/yigiter/Desktop/Dropbox/INS/SAGE/mfiles/data/';
OUTDIR='/home/yigiter/Desktop/Dropbox/INS/SAGE/mfiles/data/output/';

IMUTYPE=2;

%open files
fimu=fopen([INDIR 'imudata.bin'],'r');
imudata=fread(fimu,7,'double');
sr_a=round((TIME.interval(1)-imudata(1))/TIME.dt);
if sr_a>0
    fseek(fimu, sr_a*7*8, 0);
end

fodo=fopen([INDIR 'odo.txt'],'rt');
fgps=fopen([INDIR 'gps.txt'],'rt');

ffilres=fopen([OUTDIR 'navresult.bin'],'w');
fupdres=fopen([OUTDIR 'updresult.bin'],'w');
fpreres=fopen([OUTDIR 'preresult.bin'],'w');
fnavres=fopen([OUTDIR 'prenavresult.bin'],'w');

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

%%%Compute the alignment between the geodetic and the body frame
z_axis=-meandat(4:6)/norm(meandat(4:6));
ngyro=meandat(1:3)/norm(meandat(1:3));
vr_a=cross(z_axis,ngyro);
y_axis=vr_a/norm(vr_a);
x_axis=cross(y_axis,z_axis);
Cnb=[x_axis y_axis z_axis];
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
initime=preimutime;
covupt_time=preimutime;  %the time that we last updated the covariance
zupt_time=preimutime;   %the time that we last applied a zupt
imuaccum=zeros(6,1);
imuerrors=zeros(12,1);
gpsctr=0;

%I will accumulate 4 imu samples before processing.
imudata=fread(fimu, [7,4], 'double'); %imudata must be in increment form (not rate). If the file contains rate data, correct it here.
curimutime=imudata(1,end);

while (~feof(fimu))
    %Write result to the files   
    vr_a=posdiff_v000(Llh(1:2)', Llh(3), iniLlh);
    vr_c=quat2dcm_v000(qbn);
    fwrite(ffilres,[preimutime;vr_a';Vn;dcm2euler_v000(vr_c)*180/pi;diag(P)],'double');
    
    if (imudata(1)-initime)>60
        disp(['Process Time:' num2str(imudata(1)/60)]);
        initime=imudata(1);
    end
    
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
    
    %Check if we are going to update.
    %Hardcode the ZUPT intervals here.
    cov_upd=0;
    zupt_upd=0;
    odo_upd=0;
    gps_upd=0;
    if (curimutime-covupt_time)>=TIME.maxCovStep
        cov_upd=1;
    elseif (curimutime>ododata(1))
        cov_upd=1;
        odo_upd=1;
    elseif (curimutime>gpsdata(1))
        cov_upd=1;
        gps_upd=1;
    elseif curimutime-zupt_time>=TIME.zuptDt
        if (curimutime>TIME.zuptInterval(1,1) && curimutime<TIME.zuptInterval(1,2))
            cov_upd=1;
            zupt_upd=1;
        elseif (curimutime>TIME.zuptInterval(2,1) && curimutime<TIME.zuptInterval(2,2))
            cov_upd=1;
            zupt_upd=1;
        end
    end
    
    
    %Update the covariance
    if (cov_upd)
        %propagate the covariance
        covdt=curimutime-covupt_time;
        %[STM, Qd]=sys_ned_dcm_v001(Llh, Vn, qbn, imuaccum(4:6)/covdt, imuaccum(1:3)/covdt, covdt);
        [STM, Qd]=sys_llh_phipsi_v000(Llh, Vn, qbn, imuaccum(4:6)/covdt, imuaccum(1:3)/covdt, covdt, 1, IMUTYPE);
        P=STM*P*STM'+Qd;
        
        %Propagate the imu error states
        imuerrors=STM(10:end,10:end)*imuerrors;
        
        covupt_time=curimutime;
        imuaccum=zeros(6,1);
        
        %%Smoother Stuff
        %Record the predicted solutions
        fwrite(fnavres,[curimutime;Llh;Vn;qbn;imuerrors],'double');
        fwrite(fpreres,[curimutime;mat2vec_v000(P);STM(:);mat2vec_v000(Qd)],'double');  %I do not need Q for smoothing. But I am going to use it for disturbance estimation.
    end
    
    
    %%If there is any observation to process
    if (gps_upd || odo_upd || zupt_upd)
        %Smoother stuff
        STM_upd=eye(21);
        smo_dx_inp=0;
        smo_dP_inp=0;

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
        if (odo_upd)  %we have an odometer data
            Cnb=quat2dcm_v000(qbn)';
            if (OBS.useNonHolo)
                inno=[ododata(2);0;0]-Cnb*Vn;
                H=[zeros(3) -Cnb Cnb*skew(Vn) zeros(3,12)];
                R=diag([OBS.odoR, OBS.nonHoloR, OBS.nonHoloR]);
            else
                inno=ododata(2)-Cnb(1,:)*Vn;
                H=[zeros(1,3) -Cnb(1,:) Cnb(1,:)*skew(Vn) zeros(1,12)];
                R=OBS.odoR;
            end

            %Kalman
            klmn_scrpt;

            %Read the next odometer data
            ododata=fscanf(fodo,'%f%f',[1,2]);
            if (isempty(ododata))
                ododata(1)=inf;
            end
        end


        %%GPS updates
        if (gps_upd)
            vr_a=posdiff_v000(Llh(1:2)', Llh(3), [gpsdata(2:3)*pi/180 gpsdata(4)]);
            disp(vr_a');

            lever=diag([1/Rappx, 1/Rappx,1])*quatrot_v000(qbn,OBS.LA,0); %Thx Jianzhu
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
        mx_a=mat2vec_v000(smo_dP_inp);
        mx_b=STM_upd(:);
        fwrite(fupdres,[curimutime;smo_dx_inp;mx_a;mx_b],'double');
        
        if (gpsctr==2)
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


%Write the latest filtered result
vr_a=posdiff_v000(Llh(1:2)', Llh(3), iniLlh);
vr_c=quat2dcm_v000(qbn);
fwrite(ffilres,[preimutime;vr_a';Vn;dcm2euler_v000(vr_c)*180/pi],'double');

fclose all;
