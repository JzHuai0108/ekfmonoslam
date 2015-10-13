% this scripts is only tested on H764G and oem 615 gps data collected on
% July 19th 2012 for which it performs very well.
close all;
clear;
caseID=2;
switch caseID
    case 1
        % test on H764G 1 data
        TIME.dt=1/256;  %sampling time
        TIME.maxCovStep=1/512; %maximum covariance propagation step, smaller or equal to 1/256 means single speed mode
        OBS.LA=[ -0.746; 0.454; -1.344]; %level arm offset of gps antenna in the body and IMU frame
        
        resdir='C:\Users\huai.3\Desktop\source_testdata_specs\20130719\';
        imufile=[resdir,'H764GM0719_1.csv'];
        gpsfile=[resdir,'rtkout\oem615_20130719.pos'];
        posFileName = [resdir,'rtkout\071913_chen_Topcon2_1.pos'];
        %Time configuration
        startTime=500387.2;
        endTime=501387;
        gpsNum=300;% the number of gps data to use from the startTime
        %--------------------------------------------------------------------------
        % ZUPT start and end Time, determine before the filter
        % format n x 2
        % n is the number of segment
        % first  column is the start time
        % second column is the end time, for example [575908      576243; 576502      576602]
        options.zuptSE=[500387.2, 500600]; %Time intervals to apply zupt
        
        %NHC and ZUPT options
        rateZUPT=round(sqrt(1/TIME.dt));
        rateNHC=inf; %round(sqrt(1/TIME.dt));
        % ZUPT options
        % Sigma for ZUPT measurement constraints (<=0.0 => turn off ZUPT)
        options.sigmaZUPT = 0.1;
        
        
        % Non-Holonomic Constraint (NHC) sigma (<=0.0 => turn off NHC)
        options.sigmaNHC = 0.1;
        
        %Observation params
        OBS.useNonHolo=1;% 1;   %0:Don't use, 1:use
        OBS.zuptR=options.sigmaZUPT^2;
        
        inillh=[39.958500453*pi/180;  -83.055237599*pi/180;   205.5195];
        qbn=att2qua([-0.004486084	-0.002502441 -90.29663*pi/180]);% 500387.2 H764G INS data
        
        
    case 2
        % test on H764G 1 data collected on Aug8 2013
        TIME.dt=1/256;  %sampling time
        TIME.maxCovStep=1/512; %maximum covariance propagation step, smaller or equal to 1/256 means single speed mode
        OBS.LA=[ -0.746; 0.454; -1.344]; %level arm offset of gps antenna in the body and IMU frame
        
        
        
        isOutNED=true;
        resdir='H:\relaylatest\20130808\';
        imufile=[resdir,'H764GM0808_1.csv'];
        filresfile=[resdir, 'filresult.bin']; % navigation states
        imuresfile=[resdir, 'imuresult.bin']; % imu error terms
        kmlfilename=[resdir, 'googleearth.kml'];
        
        % imu options
        startTime=415000.0;
        % 415275 start moving,  415527.04s having camera measurements
        endTime=415400.0;
        options.imuErrorModel=4; % how bias and scale factor is modeled
        options.mechanization=2; % 1 for wander azimuth
        mode=2; % 0 for e formulation, 1 for phi, 2 for psi
        options.imutype=4;      % H764G-1
        options.dt=1/256;  %sampling interval
        options.maxCovStep=3.5*options.dt; %maximum covariance propagation step, if equal to dt, means single speed mode
        options.Cb2imu=eye(3); % set this [] in order to be estimated later
        options.Timu2body=zeros(3,1); % h764G is the body frame
        
        imuFileType=1;
        %Initial PVA of IMU
        options.inillh_ant=[40.003224712*pi/180;  -83.042989960*pi/180;   212.7058];
        inillh=options.inillh_ant;
        options.Vn=[0;0;0];
        options.Ve=[0;0;0];
        options.qb2n=att2qua([0.007873535	-0.003601074 5.394287*pi/180]);% 415000 H764G INS data
        qbn=options.qb2n;
        options.InvalidateIMUerrors=false;
        options.initAttVar=1*pi/180; % 1 deg std for roll and pitch, 2 times 5 deg for yaw std
        % gps options
        useGPS=false;
        useGPSstd=true; % use the std in the rtklib GPS solutons
        options.Tant2body=[ -0.746; 0.454; -1.344]; %level arm offset of gps antenna in the body and H764G frame
        options.gpsnum=4000*5;   % 5Hz x s, the interval of GPS coverage
        gpspostype=2;           % 1 for calender time 2 for GPSTOW format, both produced by RTKlib
        gpsfile=[resdir,'oem615_20130809.pos'];
        posFileName = gpsfile;
        % ZUPT start and end Time, determine before the filter
        % format n x 2
        % n is the number of segment
        % first  column is the start time
        % second column is the end time, for example [575908      576243; 576502      576602]
        % ZUPT options
        options.zuptSE=[415000.0, 415250.0]; %Time intervals to apply zupt
        options.sigmaZUPT = 0.1;% unit m/s
        OBS.zuptR=options.sigmaZUPT^2;
        rateZUPT=round(sqrt(1/options.dt));
        % NHC options
        options.sigmaNHC = 0.1;% unit m/s
        rateNHC=inf;% round(sqrt(1/options.dt));
        %         rateNHC=inf; % this generally cause worse results
        % minimum velocity before applying the NHC, this option decouples ZUPT and NHC
        options.minNHCVel=2.0;
        useCam=false;
        options.Cimu2cam= R2(pi/2)*R1(pi/2);
        options.Tcam2body=[2.139; -0.102; -0.925]; % casio 2 in body frame
        
        
end

% bias/sigma/data length in memory for automatically-detect ZUPT using acceleration
% measurements
options.autoZUPT_bias_fxyz = 9.8;
options.autoZUPT_sigma_fxyz = 0.01;
options.autoZUPT_bw_fxyz = 20*2;
% bias/sigma/data length in memory for automatically-detect ZUPT using gyro
% measurements
options.autoZUPT_bias_wxyz = 0.0;
options.autoZUPT_sigma_wxyz =1e-3;
options.autoZUPT_bw_wxyz = 20*2; % 20 for double sided, 20*2 for one sided forward filtering

%open files
fimu=fopen(imufile,'r');
h= fgetl(fimu);
while(true)
    if(~isempty(strfind(h,'RAW 01')))
        break;
    else h= fgetl(fimu);
    end
end
h= fgetl(fimu);
mass=textscan(h,'%f','delimiter',',');
imudata=mass{1};
while(imudata(9,1)==0||imudata(2,1)<startTime)
    preimudata=imudata(:,end);%record previous imudata
    h= fgetl(fimu);
    mass=textscan(h,'%f','delimiter',',');
    imudata=mass{1};
end
imudata=imudata([2:8],1);% gps time, xyz delta v, xyz delta theta
imudata(2:4,1)=imudata(2:4,1)*.3048;% convert to metric unit meter
preimudata=preimudata([2:8],1);% gps time, xyz delta v, xyz delta theta
preimudata(2:4,1)=preimudata(2:4,1)*.3048;% convert to metric unit meter

fgps=fopen(gpsfile,'rt');
% create the folder to store the graph and filtering results
graphdir=[resdir, 'sol\'];
if(exist(graphdir,'dir')~=7)
    mkdir(graphdir);
end

filresfile=[graphdir, 'filresult.bin'];
ffilres=fopen(filresfile,'w');
fimures=fopen([graphdir,'imuresult.bin'],'w');


%Initial PVA of IMU

Cen=llh2dcm_v000(inillh(1:2),[0,1]);
height=inillh(3); %elipsoid height
Vn=[0;0;0];


%Load IMU error definitions
IMU_ERRDEF=imu_err_defs_v000(4);

%Initial Covariance matrix
nst=21; % Rx RY RZ in earth centered n-frame, vn, RPY, acc bias, gyro bias, acc scale and gyro scale
%acc bias, gyro bias, acc scale and gyro scale are random constants in the
%filter
P=zeros(nst);
P(1:3,1:3)=eye(3)*5^2; %position errors in meter
P(4:6,4:6)=eye(3)*0.1^2; %vel error in m/s
P(7:9,7:9)=diag([1,1,3].^2*(pi/180)^2); %attitude errors in radian
P(10:12,10:12)=eye(3)*IMU_ERRDEF.acc_bias_var*100;
P(13:15,13:15)=eye(3)*IMU_ERRDEF.gyro_bias_var*100;
P(16:18,16:18)=eye(3)*IMU_ERRDEF.acc_scale_var*100;
P(19:21,19:21)=eye(3)*IMU_ERRDEF.gyro_scale_var*100;

%%%%%Discard all observation data before the current time
%remove the header of gps posdata
% h= fgetl(fgps);
% while(true)
%     if(isempty(strfind(h,'%')))
%         break;
%     else h= fgetl(fgps);
%     end
% end
% stoic=sscanf(h,'%d/%d/%d%d:%d:%f%f%f%f%d%d');
% [weeksec, weeknum]=Calender2GPSWeek(stoic(1:6));
% mess=stoic;
% gpsdata=zeros(1,7); %lattitude longitude in degree and height in meter, Q and no of satels and one reserve
% gpsdata(1)=weeksec+(mess(4:6)-stoic(4:6))'*[3600;60;1];
% gpsdata(2:6)=mess(7:11);
% while (gpsdata(1)<=preimudata(1))
%     h= fgetl(fgps);
%     mess=sscanf(h,'%d/%d/%d%d:%d:%f%f%f%f%d%d');
%     gpsdata(1)=weeksec+(mess(4:6)-stoic(4:6))'*[3600;60;1];
%     gpsdata(2:6)=mess(7:11);
% end
gpsdata=inf;
%%%%%Start the main INS
%prepare for the loop
initime=preimudata(1);
covupt_time=preimudata(1);  %the time that we last updated the covariance
imuaccum=zeros(6,1);
imuerrors=zeros(12,1);
gpsctr=0;

%256Hz.Therefore, I am going to read 1 records at a time
curimutime=imudata(1,end);
preimutime=preimudata(1,end);

%Record covariance and navigation solution for the smoother
fwrite(fimures,[curimutime;imuerrors;sqrt(diag(P(10:21,10:21)))],'double');
kkk=1;% to count how many imu epochs after the recent gps observations
while (~feof(fimu)&&curimutime<endTime)
    %Write result to the files
    vr_a=posdiff_v000(Cen(3,:), height, inillh);
    vr_c=quat2dcm_v000(qbn);
    fwrite(ffilres,[preimutime;vr_a';Vn;dcm2euler_v000(vr_c)*180/pi;sqrt(diag(P(1:9,1:9)))],'double');
    
    if (imudata(1)-initime)>60
        disp(['Process Time:' num2str(imudata(1))]);
        initime=imudata(1);
    end
    
    %Coning-sculling part (ignagni 1996/Algo7 for coning and 1998/Algo5 for
    %sculling)
    dt=curimutime-preimutime;
    gyroinc=sum(imudata(5:7,:),2);
    accinc=sum(imudata(2:4,:),2);
    angleinc=gyroinc;
    vr_a=0.5*cross(gyroinc, accinc);
    velinc=accinc+vr_a;
    
    %Calibrate the imu
    angleinc=angleinc-(imuerrors(4:6)*dt+diag(angleinc)*imuerrors(10:12)/1000);
    velinc=velinc-(imuerrors(1:3)*dt+diag(velinc)*imuerrors(7:9)/1000);
    
    %imu data accumulator for the covariance update
    imuaccum=imuaccum+[angleinc;velinc];
    
    %run strapdown with angle and velocity increments
    [qbn, Vn, Cen, height]=strapdown_Cen_quat_v000(qbn, Vn, Cen, height, velinc, angleinc, dt);
    
    %Update the covariance
    %ododata(1)=inf;
    if ((curimutime-covupt_time)>=TIME.maxCovStep) || (curimutime>gpsdata(1))
        %propagate the covariance
        covdt=curimutime-covupt_time;
        % the following two model gives very similar result though the v000
        % seems very slight little better than v002
        %         [STM, Qd]=sys_metric_phipsi_v002(Cen, height, Vn, qbn, imuaccum(4:6)/covdt, imuaccum(1:3)/covdt, covdt,2,4); %use psi implementation
        
        [STM, Qd]=sys_metric_phipsi_v000(Cen, height, Vn, qbn, imuaccum(4:6)/covdt, imuaccum(1:3)/covdt, covdt,2,4,3); %use psi implementation
        
        P=STM*P*STM'+Qd;
        
        %Propagate the imu error states
        imuerrors=STM(10:end,10:end)*imuerrors; % not necessary as we are doing close
        covupt_time=curimutime;
        imuaccum=zeros(6,1);
        
        %Record covariance and navigation solution for the smoother
        %Note that I am recording the predicted solutions. That is why the
        %backward part must use filtered solution.
        fwrite(fimures,[curimutime;imuerrors;sqrt(diag(P(10:21,10:21)))],'double');
    end
    
    isZUPT=mod(kkk, rateZUPT)==0;% isZUPT&&
    if (isZUPT&&options.sigmaZUPT>0.0&&(options.autoZUPT_bw_fxyz>=3||options.autoZUPT_bw_wxyz>=3||~isempty(options.zuptSE)))
        isStatic = false;
        %isStatic = ~isempty(find(zupt_LOC==index));
        if (~isempty(options.zuptSE) && ~isempty(find(options.zuptSE(:,1)<=curimutime && options.zuptSE(:,2)>=curimutime, 1)))
            isStatic = true;
        end
        % automatically detect ZUPT by looking forward and backward, to be
        % improved using cache imu data. I distrust automatic zupt
        % detection
        %         if (~isStatic)
        %             fxyz3 = sqrt(velinc'*velinc)/TIME.dt;
        %             rms_fxyz = abs(fxyz3-options.autoZUPT_bias_fxyz);
        %             wxyz3 = sqrt(angleinc'*angleinc)/TIME.dt;
        %             rms_wxyz = abs(wxyz3-options.autoZUPT_bias_wxyz);
        %             if (rms_wxyz<options.autoZUPT_sigma_wxyz&&rms_fxyz<options.autoZUPT_sigma_fxyz)
        %                 isStatic = true;
        %             end
        %         end
        if (isStatic)
            %--------------------------------------------------------------
            % velocity measurement update: Note Velocity in XYZ
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
    isNHC =mod(kkk, rateNHC)==0;
    if (isNHC&&options.sigmaNHC>0.0)
        %------------------------------------------------------------------
        % non-holonomic constraints
        %------------------------------------------------------------------
        alb=[0 1 0; 0 0 1];
        Cnb=quat2dcm_v000(qbn)';
        inno=alb*quatrot_v000(qbn,Vn,1);% n frame to b frame
        H=[zeros(2,3) alb*Cnb -alb*Cnb*skew(Vn) zeros(2,12)];
        R=eye(2)*options.sigmaNHC^2;
        %Kalman
        K=P*H'*inv(H*P*H'+R);
        dx=K*inno;
        P=(eye(21)-K*H)*P*(eye(21)-K*H)'+K*R*K';
        
        %Correct states
        [qbn, Vn, Cen, height]=correctnav_Cen_v000(qbn, Vn, Cen, height, dx(1:9), 2, 2);
        imuerrors=imuerrors+dx(10:end);
        %------------------------------------------------------------------
    end
    
    %%%%Apply the observations
    %%gps updates
    if (curimutime>gpsdata(1)&&(curimutime-gpsdata(1)<TIME.dt))
        %it is OK to assume spherical earth to form observations (see
        %Savage:15.1.2.2.32)
        %%%However, If you want to be a little bit more precise use the following
        kkk=0; % to count how many imu epochs after the recent gps observations
        gpsecef=ecef2geo_v000([gpsdata(2)/180*pi;gpsdata(3)/180*pi;gpsdata(4)],1);
        Re=6378137/(sqrt(1.0-0.00669437999014*Cen(3,3)^2));
        posecef=-[(Re+height)*Cen(3,1);(Re+height)*Cen(3,2);(Re*(1-0.00669437999014)+height)*Cen(3,3)];
        inno=Cen*(posecef-gpsecef)+quatrot_v000(qbn,OBS.LA,0);
        %         Rappx=6378137/(sqrt(1.0-0.00669437999014*Cen(3,3)^2))+height;
        %         gpssL=sin(gpsdata(2)*pi/180);
        %         gpscL=cos(gpsdata(2)*pi/180);
        %         gpssl=sin(gpsdata(3)*pi/180);
        %         gpscl=cos(gpsdata(3)*pi/180);
        %         vr_a=Rappx*Cen(1:2,:)*[gpscL*gpscl;gpscL*gpssl;gpssL];
        %         vr_b=quatrot_v000(qbn,OBS.LA,1);
        %         inno=[-vr_a;gpsdata(4)-height]+vr_b;
        disp(inno);
        
        H=[eye(3) zeros(3,18)]; %Completely ignore the effect of Cbn errors on level-arm compensation. If you want, use H(:,7:9) = skew(OBS.LA);
        %the following variance is ok based on the quality indicator in RTKlib output
        if(gpsdata(1,5)==1)
            R=diag([0.05,0.05,0.05].^2);
        elseif(gpsdata(1,5)==2)
            R=diag([1.0,1.0,1.0].^2);
        else
            R=diag([15,15,15].^2);
        end
        
        %Kalman
        K=P*H'*inv(H*P*H'+R);
        dx=K*inno;
        P=(eye(21)-K*H)*P*(eye(21)-K*H)'+K*R*K';
        
        %Correct states
        [qbn, Vn, Cen, height]=correctnav_Cen_v000(qbn, Vn, Cen, height, dx(1:9), 2, 2);
        imuerrors=imuerrors+dx(10:end);
        
        %Read the new gps data
        if(~feof(fgps))
            h= fgetl(fgps);
            if (~ischar(h))
                gpsdata(1)=inf;
            else
                mess=sscanf(h,'%d/%d/%d%d:%d:%f%f%f%f%d%d');
                gpsdata(1)=weeksec+(mess(4:6)-stoic(4:6))'*[3600;60;1];
                gpsdata(2:6)=mess(7:11);
            end
        else
            gpsdata(1)=inf;
        end
        gpsctr=gpsctr+1;
        if (gpsctr>=gpsNum)
            gpsdata(1)=inf;% stop using GPS, dead reckoning
        end
    end
    
    %Read the next imu data
    preimudata=imudata;
    preimutime=curimutime;
    while(preimutime>=curimutime)
        h= fgetl(fimu);
        if (~ischar(h))
            imudata=[];
            break;
        end
        mass=textscan(h,'%f','delimiter',',');
        imudata=mass{1};
        imudata=imudata([2:8],1);% gps time, xyz delta v, xyz delta theta
        imudata(2:4,1)=imudata(2:4,1)*.3048;% convert to metric unit meter
        curimutime=imudata(1,end);
    end
    if (isempty(imudata))
        break;
    else
        kkk=kkk+1;
    end
end
disp(['See you space cowboy in the backward part..']);
fclose all;
%display the results
isOutNED=true;
kf = readdata(filresfile, 1+9+9);


posdata=load(posFileName);
posdata=posdata(:,2:end);
posdata(:,2:3)=posdata(:,2:3)*pi/180;

[~, minl]=min(abs(posdata(:,1)-kf(1,1)));
if(posdata(minl,1)-kf(1,1)<0)
    minl=minl+1;
end
[~, maxl]=min(abs(posdata(:,1)-kf(end,1)));
if(posdata(maxl,1)-kf(end,1)<0)
    maxl=maxl+1;
end
maxl=min(size(posdata,1),maxl);
posdata=posdata(minl:maxl,:);
startXYZ= blh2xyz(inillh);
startCen=llh2dcm_v000(inillh(1:2),[0,1]);
display(startXYZ);
display(startCen);
for i=1:maxl-minl+1
    rovBLH = posdata(i,2:4);
    rovXYZ = blh2xyz(rovBLH);
    dXYZ = rovXYZ-startXYZ;
    dNED = startCen*dXYZ';
    posdata(i,2:4) = dNED';
end

f(21) = figure;
plot(kf(:,3),kf(:,2),'g.')
hold on
plot(posdata(:,3),posdata(:,2),'r+')
grid
axis equal
if (isOutNED)
    xlabel('East [m]')
    ylabel('North[m]')
else
    xlabel('Lontitude [radian]')
    ylabel('Latitude [radian]')
end
title('Ground Track and red GPS');
saveas(f(21),[graphdir,'red truth and track'],'fig');

f(22) = figure;

plot(kf(:,1)-kf(1,1),kf(:,4),'g.')
hold on
plot(posdata(:,1)-kf(1,1),posdata(:,4),'r+');
grid
xlabel('Time [s]')
ylabel('Height/m')
title('Height of antenna by IMU and red GPS');
saveas(f(22),[resdir 'red truth and height'],'fig')
resdir='C:\JianzhuHuai\GPS_IMU\programs\osucfm\data\20130719\';
filresfile=[graphdir, 'filresult.bin'];
kf = readdata(filresfile, 1+9+9);
plotkf(kf, graphdir);
err = readdata([graphdir, 'imuresult.bin'], 1+12+12);
ploterr(err,graphdir);
