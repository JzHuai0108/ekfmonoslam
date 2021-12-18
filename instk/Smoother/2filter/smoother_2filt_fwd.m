fclose all;
clear;

%Time configuration
TIME.interval=[0 2000];    %start/stop
TIME.dt=1/200;  %sampling time
TIME.salign=200;   %static alignment duration
TIME.zupt=[7 84];    %Time intervals to apply zupt
TIME.maxCovStep=0.1; %maximum covariance propagation step

%Observation params
OBS.odoR=0.2^2;  %(m/s)^2
OBS.useNonHolo=1;   %0:Don't use, 1:use
OBS.nonHoloR=0.1^2; %(m/s)^2
OBS.gpsR=0.1^2'; %(m)^2
OBS.zuptR=0.1^2;
OBS.LA=[0.2;0;0];

%open files
fimu=fopen('/home/yigiter/Desktop/Dropbox/INS/SAGE/mfiles/data/imudata.bin','r');
imudata=fread(fimu,7,'double');
sr_a=round((TIME.interval(1)-imudata(1))/TIME.dt);
if sr_a>0
    fseek(fimu, sr_a*7*8, 0);
end

fodo=fopen('/home/yigiter/Desktop/Dropbox/INS/SAGE/mfiles/data/odo.txt','rt');
fgps=fopen('/home/yigiter/Desktop/Dropbox/INS/SAGE/mfiles/data/gps.txt','rt');

fnavres=fopen('/home/yigiter/Desktop/Dropbox/INS/SAGE/mfiles/data/navresult.bin','w');
ffilres=fopen('/home/yigiter/Desktop/Dropbox/INS/SAGE/mfiles/data/filresult.bin','w');
fcovres=fopen('/home/yigiter/Desktop/Dropbox/INS/SAGE/mfiles/data/covresult.bin','w');

%Initial PVA
inillh=[39.88864*pi/180;32.78002*pi/180;1112];
Cen=llh2dcm_v000(inillh(1:2),[0,1]);
height=inillh(3); %elipsoid height
Vn=[0;0;0];
qbn=[1;0;0;0]; %Unknown:To be determined during coarse alignment

%Load IMU error definitions
IMU_ERRDEF=imu_err_defs_v000(1);

%Initial Covariance matrix
nst=21;
P=zeros(nst);
P(1:3,1:3)=eye(3)*5^2; %position errors in meter
P(4:6,4:6)=eye(3)*0.05^2; %vel error in m/s
P(7:9,7:9)=eye(3); %attitude errors in radian (will be determined later during coarse alignment)
P(10:12,10:12)=eye(3)*IMU_ERRDEF.acc_bias_var;
P(13:15,13:15)=eye(3)*IMU_ERRDEF.gyro_bias_var;
P(16:18,16:18)=eye(3)*IMU_ERRDEF.acc_scale_var;
P(19:21,19:21)=eye(3)*IMU_ERRDEF.gyro_scale_var;

%%%%%Static Alignment
%%Step I:Coarse align
%coarse_dur=gyro_arw/(gyro_bias_var)*2;
coarse_dur=TIME.salign;
%Compute the mean of imu data during the coarse alignment period
datbuf=fread(fimu, [7 round(coarse_dur/TIME.dt)], 'double');
meandat=mean(datbuf(2:end,:),2)/TIME.dt;
preimudata=datbuf(:,end);   %In fact, i do not need this. However, in the future I may switch to con/scull algorithm requring previous minor interval.
clear datbuf;
%Compute the alignment between the geodetic and the body frame
[Fc wen_n wie_n g]=geoparam_v001(2, Cen(:,3), height, zeros(3,1));
[Cnb E err]=align_opt_v000([meandat(4:6), meandat(1:3)],[[0;0;-g], wie_n] ,[1 0]); %Note:E defined for errors on imu. Therefore, the order of vectors is important
qbn=dcm2quat_v000(Cnb');

%%Step II:Determine i)the attitude covariance ii)cross correlations between
%%initial attitude and imu error states
Cimu=[eye(6), diag([meandat(4),meandat(5),meandat(6),meandat(1),meandat(2),meandat(3)]/1000)];
Rimu=diag([IMU_ERRDEF.acc_vrw;IMU_ERRDEF.acc_vrw;IMU_ERRDEF.acc_vrw;IMU_ERRDEF.gyro_arw;IMU_ERRDEF.gyro_arw;IMU_ERRDEF.gyro_arw]*IMU_ERRDEF.gyro_bias_var/IMU_ERRDEF.gyro_arw);
%The following is a little bit optimistic estimate (not to mention terribly
%rough.) Therefore, it is better to artifically increase Pimu(7:9,7:9) (not
%the cross correlations, just the covaraince) by arbitrarily increasing Rimu
P(7:9,7:9)=E*(Cimu*P(10:21,10:21)*Cimu'+Rimu)*E'; 
P(7:9,10:21)=E*Cimu*P(10:21,10:21);
P(10:21,7:9)=P(7:9,10:21)';

P(9,9)=(3*pi/180)^2; %%I simply do not trust the linearized error model of align_opt

%%%%%Discard all observation data before the current time
ododata=fscanf(fodo,'%f%f',[1,2]);
while (ododata(1)<=preimudata(1))
    ododata=fscanf(fodo,'%f%f',[1,2]);
end
gpsdata=fscanf(fgps,'%f%f',[1,7]);
while (gpsdata(1)<=preimudata(1))
    gpsdata=fscanf(fgps,'%f%f',[1,7]);
end

%%%%%Start the main INS
%prepare for the loop
initime=preimudata(1);
covupt_time=preimudata(1);  %the time that we last updated the covariance
imuaccum=zeros(6,1);
imuerrors=zeros(12,1);
gpsctr=0;

%I will use 4-point (single interval) coning/sculling algorithms to reduce the data rate to
%100Hz.Therefore, I am going to read 4 records at a time
imudata=fread(fimu, [7,4], 'double'); %result must be in increment form (not rate). (
curimutime=imudata(1,end);
preimutime=preimudata(1,end);

%Record covariance and navigation solution for the smoother
fwrite(fnavres,[curimutime;Cen(:);height;Vn;qbn;imuerrors],'double');
fwrite(fcovres,[curimutime;mat2vec_v000(P)],'double');

while (~feof(fimu))
    %Write result to the files   
    vr_a=posdiff_v000(Cen(3,:), height, inillh);
    vr_c=quat2dcm_v000(qbn);
    fwrite(ffilres,[preimutime;vr_a';Vn;dcm2euler_v000(vr_c)*180/pi],'double');
    
    if (imudata(1)-initime)>60
        disp(['Process Time:' num2str(imudata(1)/60)]);
        initime=imudata(1);
    end
    
    %Coning-sculling part (ignagni 1996/Algo7 for coning and 1998/Algo5 for
    %sculling)
    gyroinc=sum(imudata(2:4,:),2);
    accinc=sum(imudata(5:7,:),2);
    vr_a=cross(((54/105)*imudata(2:4,1)+(92/105)*imudata(2:4,2)+(214/105)*imudata(2:4,3)), imudata(2:4,4));
    angleinc=gyroinc+vr_a;  
    vr_a=0.5*cross(gyroinc, accinc);
    vr_b=cross(((54/105)*imudata(2:4,1)+(92/105)*imudata(2:4,2)+(214/105)*imudata(2:4,3)), imudata(5:7,4));
    vr_c=cross(((54/105)*imudata(5:7,1)+(92/105)*imudata(5:7,2)+(214/105)*imudata(5:7,3)), imudata(2:4,4));
    velinc=accinc+vr_a+vr_b+vr_c;
    dt=curimutime-preimutime;
    
    %Calibrate the imu
    angleinc=angleinc-(imuerrors(4:6)*dt+diag(angleinc)*imuerrors(10:12)/1000);
    velinc=velinc-(imuerrors(1:3)*dt+diag(velinc)*imuerrors(7:9)/1000);
    
    %imu data accumulator for the covariance update
    imuaccum=imuaccum+[angleinc;velinc];
    
    %run strapdown with angle and velocity increments
    [qbn, Vn, Cen, height]=strapdown_Cen_quat_v000(qbn, Vn, Cen, height, velinc, angleinc, dt);
    
    %Update the covariance
    %ododata(1)=inf;
    if ((curimutime-covupt_time)>=TIME.maxCovStep) || (curimutime>ododata(1))  || (curimutime>gpsdata(1))
        %propagate the covariance
        covdt=curimutime-covupt_time;
        [STM, Qd]=sys_metric_phipsi_v000(Cen, height, Vn, qbn, imuaccum(4:6)/covdt, imuaccum(1:3)/covdt, covdt,2,1); %use psi implementation
        P=STM*P*STM'+Qd;
        
        %Propagate the imu error states
        imuerrors=STM(10:end,10:end)*imuerrors;
        
        covupt_time=curimutime;
        imuaccum=zeros(6,1);
        
        %Record covariance and navigation solution for the smoother
        %Note that I am recording the predicted solutions. That is why the
        %backward part must use filtered solution.
        fwrite(fnavres,[curimutime;Cen(:);height;Vn;qbn;imuerrors],'double');
        fwrite(fcovres,[curimutime;mat2vec_v000(P)],'double');
        
        %Apply ZUPT. Don't apply zupt for each imu sample. Zupt
        %for each sample must be performed on nominal trajectory, not with
        %the Kalman filter.
        if (curimutime>12125 && curimutime<12174) || (curimutime>17250 && curimutime<17337)
            inno=Vn;
            H=[zeros(3) eye(3) zeros(3,15)];
            R=eye(3)*OBS.zuptR;
            
            %Kalman
            K=P*H'*inv(H*P*H'+R);
            dx=K*inno;
            P=(eye(21)-K*H)*P*(eye(21)-K*H)'+K*R*K';

            %Correct states
            [qbn, Vn, Cen, height]=correctnav_Cen_v000(qbn, Vn, Cen, height, dx(1:9), 2, 2);
            imuerrors=imuerrors+dx(10:end);
        end
    end
    
    %%%%Apply the observations
    %%odometer updates
    if (curimutime>ododata(1))  %we have an odometer data
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
        K=P*H'*inv(H*P*H'+R);
        dx=K*inno;
        P=(eye(21)-K*H)*P*(eye(21)-K*H)'+K*R*K';

        %Correct states
        [qbn, Vn, Cen, height]=correctnav_Cen_v000(qbn, Vn, Cen, height, dx(1:9), 2, 2);
        imuerrors=imuerrors+dx(10:end);
        
        %Read the new odometer data
        ododata=fscanf(fodo,'%f%f',[1,2]);
        if (isempty(ododata))
            ododata(1)=inf;
        end
    end
    
    %%gps updates
    if (curimutime>gpsdata(1))
        %it is OK to assume spherical earth to form observations (see
        %Savage:15.1.2.2.32)
        %%%However, If you want to be a little bit more precise use the following
        gpsecef=ecef2geo_v000([gpsdata(2)/180*pi;gpsdata(3)/180*pi;gpsdata(4)],1);
        Re=6378137/(sqrt(1.0-0.00669437999014*Cen(3,3)^2));
        posecef=-[(Re+height)*Cen(3,1);(Re+height)*Cen(3,2);(Re*(1-0.00669437999014)+height)*Cen(3,3)];
        inno=Cen*(posecef-gpsecef)+quatrot_v000(qbn,OBS.LA,0); %thx Jianzhu
%         Rappx=6378137/(sqrt(1.0-0.00669437999014*Cen(3,3)^2))+height;
%         gpssL=sin(gpsdata(2)*pi/180);
%         gpscL=cos(gpsdata(2)*pi/180);
%         gpssl=sin(gpsdata(3)*pi/180);
%         gpscl=cos(gpsdata(3)*pi/180);
%         vr_a=Rappx*Cen(1:2,:)*[gpscL*gpscl;gpscL*gpssl;gpssL];
%         vr_b=quatrot_v000(qbn,OBS.LA,0);
%         inno=[-vr_a;gpsdata(4)-height]+vr_b;
        disp(inno);
        
        H=[eye(3) zeros(3,18)]; %Completely ignore the effect of Cbn errors on level-arm compensation. If you want, use H(:,7:9) = skew(OBS.LA);
        R=eye(3)*OBS.gpsR;
        
        %Kalman
        K=P*H'*inv(H*P*H'+R);
        dx=K*inno;
        P=(eye(21)-K*H)*P*(eye(21)-K*H)'+K*R*K';

        %Correct states
        [qbn, Vn, Cen, height]=correctnav_Cen_v000(qbn, Vn, Cen, height, dx(1:9), 2, 2);
        imuerrors=imuerrors+dx(10:end);

        %Read the new gps data
        gpsdata=fscanf(fgps,'%f%f',[1,7]);
        if (isempty(gpsdata))
            gpsdata(1)=inf;
        end
        
        gpsctr=gpsctr+1;
        if (gpsctr==2)
            break;
        end
    end
    
    
    %Read the next imu data
    preimudata=imudata;
    preimutime=curimutime;
    imudata=fread(fimu, [7,4], 'double');
    if (isempty(imudata))
        break;
    end
    curimutime=imudata(1,end);
end
disp(['See you space cowboy in the backward part..']);
fclose all;
