% This program tests ekf_filter_nframe, _eframe class
% the output trajectory is of the antenna

% we face the problem that some IMU records are missing, but it is
% generally very small gap, so it is reasonable to use the received
% measurements for Kalman update once the imu epoch exceeds their epoch.

% Test cases (1) start IMU in stationary mode, do ZUPT for IMU
% for at least 1 minutes, the initial frame of IMU is s0
% frame, (2) start IMU in dynamic mode

function Test_EKF_filters()
addpath('..\instkhuai11112013'); % imu functions
addpath('..\voicebox\'); % for rotro2qr and rotqr2eu, they are more robuts
% than dcm2quat_v000 and dcm2euler_v000
import java.util.LinkedList
clear variables;
clc; close all; format longg;
fprintf('\n IMU_Kinect_Calib.m to test EKF filtering in calibrating mems IMU and Kinect!\n\n');
rng('default');

experim=2;

switch experim
    case 1
        % test on Microstrain 3dm gx3-35 data integration without camera,
        % aug 08 2013, use s0 frame integration
        % NHC is a great help for DMGX3-35 in a car, because s0 frame
        % cannot take advantage of NHC, it suffers a lot. Without NHC, s0
        % frame and n frame formulation give almost identical results.
        % use of gravity magnitude constraint does not show much improvement
        % it is advised to turn off scale factor estimate, let gravity and
        % RS2C float
        isOutNED=true;
        resdir='F:\relaylatest\20130808\';
        filresfile=[resdir, 'filresult.bin']; % navigation states
        imuresfile=[resdir, 'imuresult.bin']; % imu error terms        
        % imu options
        options.startTime=415000.0;
        options.endTime=415500.0;
        options.imuErrorModel=3; % how bias and scale factor is modeled
        options.mechanization=2; % 1 for wander azimuth     
        mode=2; % 0 for e formulation, 1 for phi, 2 for psi, 3 for local frame
        options.imutype=5; % 3dm gx3 35
        options.dt=1/100;  %sampling interval
        options.maxCovStep=options.dt; %maximum covariance propagation step, if equal to dt, means single speed mode
        options.Cb2imu=eye(3); % set this [] in order to be estimated later
        options.Timu2body=zeros(3,1); % h764G is the body frame
        options.imufile=[resdir,'microstrain_20130808imu.txt'];
         
        imuFileType=0; % 3DM GX3 -35
        %Initial PVA of IMU
        options.inillh_ant=[0.69818798263956; -1.44937359550259; 212.7058];
        % lat, log, height of the antenna, in this case, ASSUMED identical to IMU,
      
        options.imuErrors=zeros(12,1);   
        options.Vn=[0;0;0];
        options.Ve=[0;0;0];
        options.qb2n=roteu2qr('xyz',[0.007873535;	-0.003601074; 5.3888*pi/180]);
        options.InvalidateIMUerrors=false;
        options.initAttVar=1*pi/180; % 1 deg std for roll and pitch, 2 times 1 deg for yaw std
        % GPS options
        useGPS=true;
        useGPSstd=false; % use the std in the rtklib GPS solutons
        options.Tant2body=[-0.746;0.454;-1.344];
        gpsSE=[0,0]; %[options.startTime+200, options.startTime+2000]; 
        gpspostype=2;           % 1 for calender time 2 for GPSTOW format, both produced by RTKlib
        gpsfile=[resdir,'oem615_20130809.pos'];
        
        % ZUPT options
        options.sigmaZUPT = 0.05;% unit m
         % ZUPT start and end Time, determine before the filter
        % format n x 2, n is the number of segment
        % first is start time, second end time,e.g. [575908,576243; 576502,576602]
        
        zuptSE=[options.startTime, options.startTime+200];  
        rateZUPT=round(sqrt(1/options.dt)); % from tests, we see zupt has adverse effect in this case
        % NHC options
        options.sigmaNHC = 0.1;% unit m/s
        rateNHC=inf; %round(sqrt(1/options.dt));
        % minimum velocity before applying the NHC, this option decouples ZUPT and NHC
        options.minNHCVel=2.0;
        
        rateGravNorm=inf; %rateZUPT;
        sigmaGravMag=3e-3; % unit m/s^2
        
        useCam=false;
        camFile=options.imufile; % where the camera measurements come from
        options.camPoseFile=[resdir, 'kinectPose.txt']; % output camera position and attitude
        deltaPhi=[0;0;0]; % [-0.487; 1.106; -0.518]/180*pi;
        options.Cimu2cam=rotqr2ro(rvec2quat_v000(deltaPhi))*[0,0,1;-1,0,0;0,-1,0]'; % (Rc2s)'
        options.Tcam2body=[3; -6.2; -4.05]*1e-2; % unit m
        % depth camera in body frame whichStep is assumed to be the imu frame,
        
        sigmaEuler=[1; 1; 1]*pi/180; % unit rad
        sigmaTrans=[2;2;2]*1e-3; % unit m
      case 2
        % test on Microstrain 3dm gx3-35 data integration with GPS data 
        % collected by octoptor on 2014/10/08 without camera         
        % it is advised to turn off scale factor estimate
        isOutNED=true;
        resdir='C:\Users\huai.3\Desktop\huai work\OctoptorINSGPStest\temp\';
        filresfile=[resdir, 'filresult.bin']; % navigation states
        imuresfile=[resdir, 'imuresult.bin']; % imu error terms        
        % imu options
        options.startTime= 326175.9;
        options.endTime= 329732.0;
        options.imuErrorModel=4; % how bias and scale factor is modeled
        options.mechanization=2; % 1 for wander azimuth     
        mode=2; % 0 for e formulation, 1 for phi, 2 for psi, 3 for local frame s0
        options.imutype=5; % 3dm gx3 35
        options.dt=1/100;  %sampling interval
        options.maxCovStep=3*options.dt; %maximum covariance propagation step, if equal to dt, means single speed mode
        options.Cb2imu=eye(3); % set this [] in order to be estimated later
        options.Timu2body=zeros(3,1); % assume imu frame is the body frame
        options.imufile='E:\UAVdata\IMU_Microstrain\uavonvan.csv';
       
        imuFileType=4; % 3DM GX3 -35 csv
        %Initial PVA of IMU
        options.inillh_ant= [ [40.002908942;  -83.014625658]*pi/180;   190.9547];
        % lat, log, height of the antenna, in this case, ASSUMED identical to IMU,
        options.imuErrors=zeros(12,1);   
        options.Vn=[0;0;0]; % initial velocity in n frame
        options.Ve=[0;0;0]; % only used in e-frame IMU formulation
        options.Vs0=[0;0;0]; % only used in s0-frame IMU formulation
        options.qb2n= dcm2quat_v000(euler2dcm_v000([-0.02783182	-0.0120298	-0.36198288]));
        %roteu2qr('xyz',[0;	0; 72/180*pi]);
        
        options.InvalidateIMUerrors=true;
        options.initAttVar=10*pi/180; % 1 deg std for roll and pitch, 2 times 1 deg for yaw std
        % GPS options
        useGPS=true;
        useGPSstd=false; % use the std in the rtklib GPS solutons
        options.Tant2body=[8;0;118]*1e-3;
        % gps start and end time
        gpsSE=[options.startTime, options.startTime+2000];        
        gpspostype=3;           % 1 for calender time 2 for GPSTOW format, both produced by RTKlib
        gpsfile='C:\Users\huai.3\Desktop\huai work\OctoptorINSGPStest\GPSnavigation.pos';
        
        % ZUPT options
        % ZUPT start and end Time, n x 2 matrix, n is the number of segment
        % first  column is the start time, second column is the end time,
        % for example [575908      576243; 576502      576602]
        
        zuptSE=[options.startTime, options.startTime+800]; 
        options.sigmaZUPT = 0.05;% unit m/s
        rateZUPT=round(sqrt(1/options.dt)); % from tests, we see zupt has adverse effect in this case
        % NHC options
        options.sigmaNHC = 0.1;% unit m/s
        rateNHC=inf; % turn off NHC
        % minimum velocity before applying the NHC, this option decouples ZUPT and NHC
        options.minNHCVel=2.0;
        
        rateGravNorm=inf; %rateZUPT;
        sigmaGravMag=3e-3; % unit m/s^2
        
        useCam=false;
        options.useCam=useCam;
        camFile=options.imufile; % dummy here, where the camera measurements come from
        options.camPoseFile=[resdir, 'kinectPose.txt']; % output camera position and attitude
        deltaPhi=[0;0;0]; % [-0.487; 1.106; -0.518]/180*pi;
        options.Cimu2cam=rotqr2ro(rvec2quat_v000(deltaPhi))*[0,0,1;-1,0,0;0,-1,0]'; % (Rc2s)'
        options.Tcam2body=[3; -6.2; -4.05]*1e-2; % unit m
        sigmaEuler=[1; 1; 1]*pi/180; % unit rad
        sigmaTrans=[2;2;2]*1e-3; % unit m
    case 3
        % test on GPSVan ESTCP\January_21_2015 H764G-1 integration with GPS 
        % it is advised to turn off scale factor estimate
        isOutNED=true;
        resdir='C:\Users\huai.3\Desktop\huai work\OctoptorINSGPStest\temp\';
        filresfile=[resdir, 'filresult.bin']; % navigation states
        imuresfile=[resdir, 'imuresult.bin']; % imu error terms        
        % imu options
        options.startTime=329148.7178;
        options.endTime=330500; 
        options.imuErrorModel=4; % how bias and scale factor is modeled
        options.mechanization=2; % 1 for wander azimuth     
        mode=1; % 0 for e formulation, 1 for phi, 2 for psi, 3 for local frame
        options.imutype=4; % h764g
        options.dt=1/256;  %sampling interval
        options.maxCovStep=3*options.dt; %maximum covariance propagation step, if equal to dt, means single speed mode
        options.Cb2imu=eye(3); % set this [] in order to be estimated later
        options.Timu2body=zeros(3,1); % assume imu frame is the body frame
        options.imufile='F:\GM0121_1.csv';       
        imuFileType=1; % h764g csv
        %Initial PVA of IMU
        options.inillh_ant=ecef2geo_v000([592601.6449  -4856620.7550   4078398.4214]',0);
        % lat, log, height of the antenna, in this case, ASSUMED identical to IMU
        options.imuErrors=zeros(12,1);   
        options.Vn=[0;0;0]; % initial velocity in n frame
        options.Ve=[0;0;0]; % only used in e-frame IMU formulation
        options.Vs0=[0;0;0]; % only used in s0-frame IMU formulation
        options.qb2n=roteu2qr('xyz',[0.00769043;	0.01144409;	-83.95203*pi/180]);

        options.InvalidateIMUerrors=true;
        options.initAttVar=2*pi/180; % 1 deg std for roll and pitch, 2 times 1 deg for yaw std
        % GPS options
        useGPS=true;
        useGPSstd=true; % use the std in the rtklib GPS solutons
        options.Tant2body=[ -0.746+1.92; 0.454; -1.344-0.5];
        % gps start and end time
        gpsSE=[options.startTime+700, options.startTime+2000];        
        gpspostype=3;           % 1 for calender time 2 for GPSTOW format, both produced by RTKlib
        gpsfile='F:\OEM615-van-front.pos';
        
        % ZUPT options
        % ZUPT start and end Time, n x 2 matrix, n is the number of segment
        % first  column is the start time, second column is the end time,
        % for example [575908      576243; 576502      576602]
        
       zuptSE = [0 329886.5; 329902.5 330059;
...

330087 330114; 330185 330209; 330282 330307;
...

330704 330774.5; 330969 331001; 331052 331066;
...

331075 331299; 331314 331449; 331461 331570.5;
...

331595 331706; 331716 331812; 331816 331831.5;
...

331844 331956.5; 331998 332108.5; 332121 332234;
...

332244 332376; 332405 332609; 333066 333278];
        
        options.sigmaZUPT = 0.1;% unit m/s
        rateZUPT=round(sqrt(1/options.dt)); % from tests, we see zupt has adverse effect in this case
        % NHC options
        options.sigmaNHC = 0.1;% unit m/s
        rateNHC=inf; % turn off NHC
        % minimum velocity before applying the NHC, this option decouples ZUPT and NHC
        options.minNHCVel=2.0;
        
        useCam=false;      % the following is dummy section for camera
        options.useCam=useCam;
        camFile=options.imufile; % dummy here, where the camera measurements come from
        options.camPoseFile=[resdir, 'kinectPose.txt']; % output camera position and attitude
        deltaPhi=[0;0;0]; % [-0.487; 1.106; -0.518]/180*pi;
        options.Cimu2cam=rotqr2ro(rvec2quat_v000(deltaPhi))*[0,0,1;-1,0,0;0,-1,0]'; % (Rc2s)'
        options.Tcam2body=[3; -6.2; -4.05]*1e-2; % unit m
        sigmaEuler=[1; 1; 1]*pi/180; % unit rad
        sigmaTrans=[2;2;2]*1e-3; % unit m
    otherwise
        error('Unsupported testing case!');
end
Counter.numimurecords=40; % the number of data put in preimudata
% Initialize the model state and covariance of state, process noise and
% measurment noise
switch(mode)
    case 0
        filter =EKF_filter_eframe(options);
    case {1,2}
        filter =EKF_filter_nframe(options);
    case 3
        filter =EKF_filter_s0framelet(options);
end

preimudata=LinkedList();% record the previous imu data

% read in imu data
[fimu, imudata, preimudata]=readimuheader(options.imufile, preimudata, options.startTime, Counter, imuFileType);
lastimu=preimudata.getLast();
preimutime=lastimu(1,end);
imuctr=1;   % to count how many IMU data after the latest GPS observations
% read the GPS data and align the GPS data with the imu data
if(useGPS)
    [fgps, gpsdata]=readgpsheader(gpsfile, gpsSE(1,1), gpspostype);
else gpsdata=inf;
end
gpsctr=1;       % the number of GPS data that has been used

%load the camera measurements with timestamp
if(useCam)
    camdata=load_KinectDataset(camFile, options.startTime);
    euc2c0imu=zeros(size(camdata,1),4);
    euc2c0imu(:,1)=camdata(:,1);
    imgepoch=camdata(1,1);
    whichStep=1;
else imgepoch=inf;
end
% Start the main INS
initime=preimutime;
imuaccum=zeros(6,1);    % record the accumulated imu measurements
curimutime=imudata(1,end);
covupt_time=preimutime; %the time that we last updated the covariance

fprintf('Estimating trajectory...\n');
ffilres=fopen(filresfile,'Wb'); % navigation states, 'W' use buffer and binary fwrite()
fimures=fopen(imuresfile,'Wb'); % imu errors

while (~feof(fimu)&&curimutime<options.endTime)
    %Write IMU's position and velocity, rpy and accel and gyro bias to the files
    if(isOutNED)
        filter.SaveToFile(options.inillh_ant, preimutime, ffilres);
    else
        filter.SaveToFile([], preimutime, ffilres);
    end
    if (imudata(1)-initime)>60
        disp(['Process Time:' num2str(imudata(1))]);
        initime=imudata(1);
    end
    % time propagation of the IMU
    imuaccum=filter.ffun_state(imuaccum,[imudata(2:7);preimutime;curimutime]);
    
    %Update the covariance
    if ((curimutime-covupt_time)>=options.maxCovStep)
        %propagate the covariance
        filter.ffun_covariance(imuaccum, covupt_time, curimutime );
        covupt_time=curimutime;
        imuaccum=zeros(6,1);
        %Record covariance and navigation solution for the smoother
        %Note that I am recording the predicted solutions. That is why the
        %backward part must use filtered solution.
        imuerrors=filter.imuErrors;
        adjcoeff=ones(size(imuerrors));
        % since the scale factor is in ppt rather than in ppm, we need to
        % convert it to ppm
        adjcoeff(filter.imuScaleFactorSIP-filter.imuBiasDriftSIP+(1:6))=1e-3;
        outimuerrors=imuerrors.*adjcoeff;
        fwrite(fimures,[curimutime;outimuerrors;...
            full(sqrt(diag(filter.p_k_k(filter.imuBiasDriftSIP+(0:11),...
            filter.imuBiasDriftSIP+(0:11))))).*adjcoeff],'double');
    end
    %Apply ZUPT. Don't apply zupt for each imu sample. Zupt
    %for each sample must be performed on nominal trajectory, not with
    %the Kalman filter.
    isStatic =~isempty(zuptSE) && ~isempty(find(((zuptSE(:,1)<=curimutime)&(zuptSE(:,2)>=curimutime))==1,1));
    isZUPT =mod(imuctr, rateZUPT)==0;
    if (isStatic&&isZUPT)
        measure=zeros(3,1);
        if(strcmp(filter.tag, 'EKF_IMU_GPS_EFRM'))
            predict=filter.rvqs2e(4:6);
            H=[zeros(3) eye(3) zeros(3,size(filter.p_k_k,1)-6)];
        elseif(strcmp(filter.tag, 'EKF_IMU_GPS_NFRM'))
            predict=filter.Vn;
            H=[zeros(3) eye(3) zeros(3, filter.covDim-6)];
        else        
        predict=filter.rvqs0(4:6);
        H=sparse([zeros(3) eye(3) zeros(3,size(filter.p_k_k,1)-6)]);
        end
        R=eye(3)*options.sigmaZUPT^2;
        filter.correctstates(predict,measure, H,R);
    end
    % apply gravity norm measurement
%     isGravNorm =mod(imuctr, rateGravNorm)==0;
%     if (isStatic&&isGravNorm)
%         %Gravity (most time consuming part of ecef implementations)
%         xyz_imu=filter.rqs02e(1:3)+quatrot_v000(filter.rqs02e(4:7),filter.rvqs0(1:3),0);
%         Llh=ecef2geo_v000(xyz_imu,0);
%         [Rn, Re, gn, sL, cL, WIE_E]=geoparam_v000(Llh);
%         predict=sqrt(imudata(2:4)'*imudata(2:4));
%         gs0=quatrot_v000(filter.rvqs0(7:10), imudata(2:4),1); % predicted gravity in s0 frame
%         H=[zeros(1,9) 1/predict*gs0' zeros(1,filter.covDim-12)];
%         measure=norm(gn,2);
%         R=sigmaGravMag^2;
%         filter.correctstates(predict,measure, H,R);
%     end
    isNHC =mod(imuctr, rateNHC)==0;
    velnorm=filter.GetVelocityMag();
    if (isNHC&&velnorm>options.minNHCVel)
        % non-holonomic constraints
        if(strcmp(filter.tag, 'EKF_IMU_GPS_EFRM'))
            Cs2e=quat2dcm_v000(filter.rvqs2e(7:10));
            Ce2b=(Cs2e*filter.Cb2imu)';
            curvel=filter.rvqs2e(4:6);
            predict=Ce2b(2:3,:)*curvel;
            H=[zeros(2,3) Ce2b(2:3,:) -Ce2b(2:3,:)*skew(curvel) zeros(2,size(filter.p_k_k,1)-9)];
        elseif(strcmp(filter.tag, 'EKF_IMU_GPS_NFRM'))
            Cs2n=quat2dcm_v000(filter.qs2n);
            Cn2b=(Cs2n*filter.Cb2imu)';
            curvel=filter.Vn;
            predict=Cn2b(2:3,:)*curvel;
            H=[zeros(2,3) Cn2b(2:3,:) -Cn2b(2:3,:)*skew(curvel) zeros(2,filter.covDim-9)];
        else            
            Cs02s=quat2dcm_v000(filter.rvqs0(7:10));
            Cs02b=filter.Cb2imu'*Cs02s;
            predict=Cs02b(2:3,:)*filter.rvqs0(4:6);
            larry=filter.Cb2imu'*skew(quatrot_v000(filter.rvqs0(7:10),filter.rvqs0(4:6),0));
            H=sparse([zeros(2,9) Cs02b(2:3,:) larry(2:3,:) zeros(2,size(filter.p_k_k,1)-15)]);
        end
        measure=[0;0];
        R=eye(2)*options.sigmaNHC^2;
        filter.correctstates(predict,measure, H,R);
    end
 
    %%Apply the  GPS observations    
    gpsSErow=find(((gpsSE(:,1)<=gpsdata(1))&(gpsSE(:,2)>=gpsdata(1)))==1,1);
    if (~isempty(gpsSErow) && curimutime>=gpsdata(1))
        imuctr=0; % to count how many imu epochs after the recent gps observations
        if(gpspostype<3)
            gpsecef=ecef2geo_v000([gpsdata(2)/180*pi;gpsdata(3)/180*pi;gpsdata(4)],1);
        else
            gpsecef=gpsdata(2:4);
        end
        measure=gpsecef;        
        if(strcmp(filter.tag, 'EKF_IMU_GPS_EFRM'))
            lever=quatrot_v000(filter.rvqs2e(7:10),filter.Tant2imu,0);
            predict=filter.rvqs2e(1:3)+lever;
            H=[eye(3) zeros(3) skew(lever) zeros(3,size(filter.p_k_k,1)-9)];
        elseif(strcmp(filter.tag, 'EKF_IMU_GPS_NFRM'))            
            Ce2n=filter.Cen;
            heights2n=filter.height;
            measure=Ce2n*gpsecef;            
            Re=6378137/(sqrt(1.0-0.00669437999014*Ce2n(3,3)^2));
            posecef=-[(Re+heights2n)*Ce2n(3,1);(Re+heights2n)*Ce2n(3,2);(Re*(1-0.00669437999014)+heights2n)*Ce2n(3,3)];
            lever=quatrot_v000(filter.qs2n,filter.Tant2imu,0);
            predict=Ce2n*posecef+lever;
            H=[eye(3) zeros(3) skew(lever) zeros(3,filter.covDim-9)];
        else            
            lever=quatrot_v000(filter.rvqs2e(7:10),filter.Tant2imu,0);
            predict=filter.rvqs2e(1:3)+lever;
            H=sparse([eye(3), skew(quatrot_v000(filter.rqs02e(4:7),filter.rvqs0(1:3),0)+lever),...
                quat2dcm_v000(filter.rqs02e(4:7)),zeros(3), ...
                -quat2dcm_v000(filter.rvqs2e(7:10))*skew(filter.Tant2imu), zeros(3,size(filter.p_k_k,1)-15)]);
        end
        % the following setting of noise variances is suitable for RTKlib output
        if(~useGPSstd)
            if(gpsdata(5)==1)
                R=diag([0.05,0.05,0.15].^2);
            elseif(gpsdata(5)==2)
                R=diag([1.0,1.0,2.0].^2);
            else
                R=diag([15,15,15].^2);
            end
            if(~strcmp(filter.tag, 'EKF_IMU_GPS_NFRM'))
                Ce2n=llh2dcm_v000(ecef2geo_v000(filter.rvqs2e(1:3,1),0),[0;1]);
                R=Ce2n'*R*Ce2n;
            end
        else
            R=4*diag(gpsdata(7:9).^2);
            if(gpspostype<3&&(~strcmp(filter.tag, 'EKF_IMU_GPS_NFRM')))
                Ce2n=llh2dcm_v000(ecef2geo_v000(filter.rvqs2e(1:3,1),0),[0;1]);
                R=Ce2n'*R*Ce2n;
            elseif(gpspostype>2&&(strcmp(filter.tag, 'EKF_IMU_GPS_NFRM')))
                R=filter.Cen*R*filter.Cen';
            end
        end
        filter.correctstates(predict,measure, H,R);        
        %Read the next gps data
        [fgps, gpsdata]=grabnextgpsdata(fgps, gpspostype);
       
        if (gpsdata(1)>gpsSE(gpsSErow,2))
            disp(['GPS outage starts from ' num2str(gpsdata(1)) 'GTOW sec!']);
            gpsdata(1)=inf;% stop using GPS, dead reckoning
        end
    end
    % image measurement
    if(curimutime>=imgepoch)
        % there are some gap in the imu recordings
        if(abs(curimutime-imgepoch)>options.dt)
            disp(['abs(curimutime-imgepoch)>options.dt at ' num2str(curimutime)...
                ' GTOW sec and image index ' num2str(whichStep) '!']);
        end
        % predict the image cooridnates of points in the states
        RTcnc0=camdata(whichStep,2:7)';
        % 6x1 vector. RTcnc0(1:3) euler angle of Rc2c0, RTcnc0(4:6), Tc0 2c
        [Rc2c0, Tc02c , H]=filter.computeH();
        predict=[unskew(eye(3)-Rc2c0*(roteu2ro('xyz',RTcnc0(1:3)))', 0); Tc02c-RTcnc0(4:6)];
        euc2c0imu(whichStep, :)=[imgepoch; rotro2eu('xyz', Rc2c0)]';
        % predicted error =predicted value-measured value,
        % measured error= measured value- measured value= 0, dummy
        measure=zeros(6,1); % dummy
        
        R=diag([sigmaEuler.^2;sigmaTrans.^2]);
        filter.correctstates(predict,measure, H,R);
        filter.SaveCamPoseandRqs02e(imgepoch,whichStep);
        
        % take the next image time
        if(whichStep<size(camdata,1))
            whichStep=whichStep+1;
            imgepoch=camdata(whichStep,1);
        else
            imgepoch=inf; % stop using images
        end
    end
    %     filter.SaveCamPoseandRqs02e(curimutime,whichStep);
    %Read the next imu data
    if(preimudata.size()==Counter.numimurecords)
        preimudata.removeFirst();
    end
    preimudata.addLast(imudata(:,end));    % record previous imudata
    preimutime=curimutime;
    [fimu, imudata]=grabnextimudata(fimu, preimutime, imuFileType);
    if (isempty(imudata))
        break;
    else
        curimutime=imudata(1,end);
        imuctr=imuctr+1;
    end
end
fprintf(' done.\n\n');
fclose all;
% display the navigation results
plotnav;
end