% This program implements the monocular IMU GPS integration based on Jones and Soatto, 2010.
% Formulated in global frame without camera input, and in local frame when
% camera input comes. KLT tracker is used for point tracking
% This program uses three classes, KLTtracker class, the EKFfilter class
% and the loopCloser class, all classes are tailored for this application.
% care should be taken when reuse them.

% data structure member in LoopCloser class
% feature base stores all manifest feature appearances in the memory

% the output trajectory is that of the antenna for it avoids transforming
% the antenna position in the rtklib solution to imu position in the plots
% of comparison section

% we face the problem that some IMU records are missing, but it is
% generally very small gap, so it is reasonable to use the received
% measurements for Kalman update once the imu epoch exceeds their epoch.

% Test cases (1) verify given 0 video frames (useCam=false) and 10 video frames to
% the filter do not make a difference in the outcomes,
% (2) two speed state update and covariance update, set model.maxCovStep>model.dt
% (3) static mode, random noise as imu input, run the imu and image processing
% components
% (4) globally initialized, without GPS, run CAM+IMU
% (5) run CAM+IMU in static situation
% (6) normal imu input and normal image input, the test image sequences of
% SPIN lab are used to verify EKF SLAM
% (7) output antenna local NED trajectory or antenna global llh trajectory, set
% isOutNED=false
% (8) no gps data, only initial values for static H764G IMU does not
% converge unless ZUPT is applied. In this almost free inertial mode, NHC
% can have a little adverse effect.
addpath('..\instk'); % imu functions
addpath('..\utilities'); % imu functions
addpath('..\toolbox_calib');% to use project_point2.m and normalize.m
addpath('..\mexopencv\');% cvcalcopticflowpyrlk and cvgoodfeaturestotrack
addpath('..\voicebox\'); % for rotro2qr and rotqr2eu, they are more robust
% than dcm2quat_v000 and dcm2euler_v000
addpath('..\EKF_monoSLAM_1pRANSAC\matlab_code\'); % ekf slam camera functions of civera

import java.util.LinkedList
clear variables;
clc; close all; format longg;
fprintf('\n mono_slam.m to test EKF filtering in combining GPS, mems IMU and monocular!\n\n');
rng('default');

experim=0; %test case for GPS IMU CAMera intergation

switch experim
    case 0
        % aug 08 2013 data MEMS 3DM GX 3 35 with bicycling people on the
        % scene, using gps solution by rtklib by Zhang Xi
        % The model we used in this case considers combination of GPS, IMU,
        % CAM, ZUPT and NHC
        % in this experiment, the MEMS 3dm gx 3-35 is mounted rigidly on
        % the vehicle. with coarse alignment, it is estimated that
        % Cimu2body has rotation angles less than 1 degree on all 3 axes.
        % output settings
        isOutNED=true;
        datadir='.\data\20130808\';
        resdir='.\data\20130808\temp\';
        
        filresfile=[resdir, 'filresult.bin']; % navigation states
        imuresfile=[resdir, 'imuresult.bin']; % imu error terms
        kmlfilename=[resdir, 'googleearth.kml'];
        
        % the input options
        % IMU options
        options.startTime=415500.00; 
        % 415275 start moving,  415527.04s having camera measurements
        options.endTime=415860.00; % 5 min after GPS outage
        options.imuErrorModel=3; % how bias and scale factor is modeled
        options.mechanization=2; % 1 for wander azimuth
        mode= 0; % 0 for e formulation, 1 for phi, 2 for psi, 3 for s0 local frame
        % it is observed that e formulation is simple to write but n frame formulation is robust
        options.imutype=5;      % MEMS 3DM GX 3-35
        options.dt=1/100;
        % maximum covariance propagation step, this is
        % set to determine the interval for covariance update
        % it is observed that maxCovStep=dt~5xdt give visually identical results
        options.maxCovStep=1/20;
        options.Cb2imu=eye(3);
        options.Timu2body=[0.03;0.03; -0.06]; % 3dm gx 3-35 in body frame
        % imufile may also be used for coarse alignment to estimate Cb2n
        options.imufile=[datadir 'microstrain_20130808imu.txt']; 
        imuFileType=0;
        options.imuErrors=zeros(12,1);
        % Initial PVA, note this is the GPS antenna position, not IMU
        %         options.inillh_ant=[40.003224712*pi/180;  -83.042989960*pi/180  ;212.7058];%415000
        options.inillh_ant=[40.003784368*pi/180  -83.042904239*pi/180   212.5730 ]'; % 415500
        options.Vn=[0.492035320491468;-4.93569659962911;-0.102610975136507];%  H764G+GPS reference at 415500
%         options.Ve=[0;0;0]; % velocity of the MEMS chassis or sensor frame in e frame
        options.Ve=llh2dcm_v000(options.inillh_ant,[0;1])'*options.Vn; 
        %         options.qb2n=att2qua([0.007965088	-0.003631592 5.39978]/180*pi);% H764G INS data, 414840
        options.qb2n=att2qua([0.752710886563209; 1.11608842515955; -84.8149334249363]/180*pi);% Roll pitch, AZ, H764G+GPS reference at 415500
        options.InvalidateIMUerrors=true;
        options.initAttVar=3*pi/180; % 5 deg std for roll and pitch, 2 times 5 deg for yaw std       
        
        % GPS options
        useGPS=true;
        useGPSstd=false; % use the std in the rtklib GPS solutons
        options.Tant2body=[-0.746;0.454;-1.344];
        options.gpsnum=100*5;   % 5Hz x s, the interval of GPS coverage
        gpspostype=2;           % 1 for calender time 2 for GPSTOW format, both produced by RTKlib
        gpsfile=[datadir,'oem615_20130809.pos'];
     
        % ZUPT options
        options.zuptSE=[414840, 415275-100];% zupt start and end time
        options.sigmaZUPT = 0.1;% unit m/s
        rateZUPT=round(sqrt(1/options.dt));
        
        % NHC options
        options.sigmaNHC = 0.3;% unit m/s
        rateNHC=round(sqrt(1/options.dt));
%         rateNHC=inf; % this generally cause worse results
        % minimum velocity before applying the NHC, this option decouples ZUPT and NHC
        options.minNHCVel=2.0;
        
        % Camera options
        options.useCam=true;
        options.camtype=1;                  % 2 Nikon D800, 1 for casio 2
        options.Cimu2cam= R2(pi/2)*R1(pi/2);
        % R1(pi/2)*R3(pi/2);% the mems frame to the camera frame
        % this equals to roteu2ro('xy', [-pi/2; -pi/2]) and R2(pi/2)*R1(pi/2)
        options.Tcam2body=[2.139; -0.102; -0.925]; % casio 2 in body frame
        
        options.sigmaCAM =1;                % unit pixel
        maxEdge=1000;% resize to less than maxEdge for each edge
        ftTblFile=[resdir, 'featureTable.txt']; % output feature table of the tracker
        options.camPoseFile=[resdir, 'camPose.txt']; % output q s2c and Ts in c, camera calibration parameters
        % this file generated by maximizing the correlation between the
        % rotation angle increments obtained from camera and imu recordings
        imgtimefile=[datadir, 'casio2timestamps.txt'];
        sequencePath =[datadir, 'casio2_3430.MOV'];
        videoName=[resdir, 'trackres.avi']; % output video
        imgseqtype=1; % 0 for image sequence and 1 for video
        initIm =12180; % 415530 GTOW
        lastIm =22000; % at least 5 minutes
   
    case 2
        % test on H764G 1 data integration with GPS, July 19 2013
        % output options
        isOutNED=true;
        resdir='C:\Users\huai.3\Desktop\source_testdata_specs\20130719\';
        filresfile=[resdir, 'filresult.bin']; % navigation states
        imuresfile=[resdir, 'imuresult.bin']; % imu error terms
        kmlfilename=[resdir, 'googleearth.kml'];
        
        % imu options
        options.startTime=500387.2;
        options.endTime=500887;
        options.imuErrorModel=4; % how bias and scale factor is modeled, random constant usually better than random walk
        options.mechanization=2; % 1 for wander azimuth
        mode=2; % 0 for e formulation, 1 for phi, 2 for psi
        options.imutype=4;      % H764G-1
        options.dt=1/256;  %sampling time
        options.maxCovStep=3.5*options.dt/2; %maximum covariance propagation step, if equal to dt, means single speed mode
        options.Cb2imu=eye(3);
        options.Timu2body=zeros(3,1); % h764G is the body frame
        options.imufile=[resdir,'H764GM0719_1.csv'];
        imuFileType=1;
        options.imuErrors=zeros(12,1);
        %Initial PVA of IMU
        options.inillh_ant=[39.958500453*pi/180;  -83.055237599*pi/180;   205.5195];
        options.Vn=[0;0;0];
        options.Ve=[0;0;0];
        options.qb2n=att2qua([-0.004486084	-0.002502441 -90.29663*pi/180]);% 500387.2 H764G INS data
        options.InvalidateIMUerrors=false; % for H764G false is doing well
        options.initAttVar=2*pi/180; % 5 deg std for roll and pitch, 2 times 5 deg for yaw std
        % gps options
        useGPS=true;
        useGPSstd=true; % use the std in the rtklib GPS solutons
        options.Tant2body=[ -0.746; 0.454; -1.344]; %level arm offset of gps antenna in the body and H764G frame
        options.gpsnum=10*5;   % 5Hz x s, the interval of GPS coverage
        gpspostype=1;           % 1 for calender time 2 for GPSTOW format, both produced by RTKlib
        gpsfile=[resdir,'rtkout\oem615_20130719.pos'];
    
        % ZUPT start and end Time, determine before the filter
        % format n x 2
        % n is the number of segment
        % first  column is the start time
        % second column is the end time, for example [575908      576243; 576502      576602]
        % ZUPT options
        options.zuptSE=[500387.2, 500600]; %Time intervals to apply zupt, 
        % if we do not use any measurements, then with good initialization,
        % the filter drifts off, but if we only use ZUPT, then the solution
        % agrees with the spreadsheet. This verifies the importance of ZUPT
        options.sigmaZUPT = 0.1;% unit m/s
        rateZUPT=round(sqrt(1/options.dt));
        % NHC options
        options.sigmaNHC = 0.1;% unit m/s
        rateNHC=inf; %round(sqrt(1/options.dt));
        %         rateNHC=inf; % this generally cause worse results
        % minimum velocity before applying the NHC, this option decouples ZUPT and NHC
        options.minNHCVel=2;
        options.useCam=false;
        options.Cimu2cam= R2(pi/2)*R1(pi/2);
        options.Tcam2body=[2.139; -0.102; -0.925]; % casio 2 in body frame
    case 3
        % test on H764G 1 data integration with GPS, Aug 8 2013
        % output options
        isOutNED=false;
        resdir='F:\relaylatest\20130808\';
        filresfile=[resdir, 'filresult.bin']; % navigation states
        imuresfile=[resdir, 'imuresult.bin']; % imu error terms
        kmlfilename=[resdir, 'googleearth.kml'];
        
        % imu options
        options.startTime=415000.0;
        % 415275 start moving,  415527.04s having camera measurements
        options.endTime=416000.0;    
        options.imuErrorModel=3; % how bias and scale factor is modeled
        options.mechanization=2; % 1 for wander azimuth
        mode=1; % 0 for e formulation, 1 for phi, 2 for psi
        options.imutype=4;      % H764G-1
        options.dt=1/256;  %sampling interval
        options.maxCovStep=3.5*options.dt; %maximum covariance propagation step, if equal to dt, means single speed mode
        options.Cb2imu=eye(3);
        options.Timu2body=zeros(3,1); % h764G is the body frame
        options.imufile=[resdir,'H764GM0808_1.csv'];
        imuFileType=1;
        options.imuErrors=zeros(12,1);
        %Initial PVA of IMU
        options.inillh_ant=[40.003224712*pi/180;  -83.042989960*pi/180;   212.7058];
        options.Vn=[0;0;0];
        options.Ve=[0;0;0];
        options.qb2n=att2qua([0.007873535	-0.003601074 5.394287*pi/180]);% 415000 H764G INS data
        options.InvalidateIMUerrors=false;
        options.initAttVar=1*pi/180; % 1 deg std for roll and pitch, 2 times 5 deg for yaw std
        % gps options
        useGPS=true;
        useGPSstd=true; % use the std in the rtklib GPS solutons
        options.Tant2body=[ -0.746; 0.454; -1.344]; %level arm offset of gps antenna in the body and H764G frame
        options.gpsnum=40*5;   % 5Hz x s, the interval of GPS coverage
        gpspostype=2;           % 1 for calender time 2 for GPSTOW format, both produced by RTKlib
        gpsfile=[resdir,'oem615_20130809.pos'];
    
        % ZUPT start and end Time, determine before the filter
        % format n x 2
        % n is the number of segment
        % first  column is the start time
        % second column is the end time, for example [575908      576243; 576502      576602]
        % ZUPT options
        options.zuptSE=[415000.0, 415250.0]; %Time intervals to apply zupt
        options.sigmaZUPT = 0.1;% unit m/s
        rateZUPT=round(sqrt(1/options.dt));
        % NHC options
        options.sigmaNHC = 0.1;% unit m/s
        rateNHC=round(sqrt(1/options.dt));
        %         rateNHC=inf; % this generally cause worse results
        % minimum velocity before applying the NHC, this option decouples ZUPT and NHC
        options.minNHCVel=2.0;
        options.useCam=false;
        options.Cimu2cam= R2(pi/2)*R1(pi/2);
        options.Tcam2body=[2.139; -0.102; -0.925]; % casio 2 in body frame
        
    case 4
        % test on Steval MKI062V2 data integration with camera, Mar 2014
        % output options
        isOutNED=true;
        resdir='F:\OpenShoe-Matlab-Implementation\';
        filresfile=[resdir, 'filresult.bin']; % navigation states
        imuresfile=[resdir, 'imuresult.bin']; % imu error terms
        kmlfilename=[resdir, 'googleearth.kml'];
        
        % imu options
        options.startTime=1310.76;    
        options.endTime=1380.0;    
        options.imuErrorModel=3; % how bias and scale factor is modeled
        options.mechanization=2; % 1 for wander azimuth
        mode=2; % 0 for e formulation, 1 for phi, 2 for psi
        options.imutype=6;      %  Steval MKI062V2 
        options.dt=1/50;  %sampling interval
        options.maxCovStep=options.dt; %maximum covariance propagation step, if equal to dt, means single speed mode
        options.Cb2imu=eye(3);
        options.Timu2body=zeros(3,1); % h764G is the body frame
        options.imufile=[resdir,'PoseIMUKinect3.tsv'];
        imuFileType=2;
        options.imuErrors=zeros(12,1);
        options.imuErrors(4:6)=-[ -0.0084376;       -0.057877;       -0.068204]; % for the Steval case, comment this line if not steval
        options.imuErrors(1:3)=[ -0.2198;   -0.2310;    0.9075];  % for the Steval case, comment this line if not steval
        
        
        %Initial PVA of IMU
        options.inillh_ant=[40.00311241687*pi/180; -83.01529294250*pi/180; 227.901]; 
        % lat, log, height guessed from http://www.daftlogic.com/sandbox-google-maps-find-altitude.htm
      
        options.Vn=[0;0;0];
        options.Ve=[0;0;0];
        options.qb2n=rotro2qr(R3(0.1477)*[0, -1,0; -1,0, 0; 0,0, -1]);
        options.InvalidateIMUerrors=false;
        options.initAttVar=3*pi/180; % 1 deg std for roll and pitch, 2 times 5 deg for yaw std
        % gps options
        useGPS=false;
        useGPSstd=true; % use the std in the rtklib GPS solutons
        options.Tant2body=[ -0.746; 0.454; -1.344]; %level arm offset of gps antenna in the body and H764G frame
        options.gpsnum=4000*5;   % 5Hz x s, the interval of GPS coverage
        gpspostype=2;           % 1 for calender time 2 for GPSTOW format, both produced by RTKlib
        gpsfile=[resdir,'oem615_20130809.pos'];
    
        % ZUPT start and end Time, determine before the filter
        % format n x 2
        % n is the number of segment
        % first  column is the start time
        % second column is the end time, for example [575908      576243; 576502      576602]
        % ZUPT options
        options.zuptSE=[1310.72, 1335.00]; %Time intervals to apply zupt
        options.sigmaZUPT = 0.05;% unit m/s
        rateZUPT=round(sqrt(1/options.dt));
        % NHC options
        options.sigmaNHC = 0.1;% unit m/s
        rateNHC=inf; %round(sqrt(1/options.dt));
        %         rateNHC=inf; % this generally cause worse results
        % minimum velocity before applying the NHC, this option decouples ZUPT and NHC
        options.minNHCVel=2.0;
        options.useCam=false;
        options.Cimu2cam= R2(pi/2)*R1(pi/2);
        options.Tcam2body=[2.139; -0.102; -0.925]; % casio 2 in body frame
    otherwise
        error('Unsupported testing case!');
end

if exist(resdir, 'dir') ~= 7
     mkdir(resdir);
end

Counter.numimurecords=40; % the number of data put in preimudata
Counter.numcamconfigrecords=300; % number of camera configuration records in camconfighistory
% Initialize the model state and covariance of state, process noise and
% measurment noise
switch(mode)
    case 0
        filter =EKF_filter_eframe(options);
    case {1,2}
        filter =EKF_filter_nframe(options);
end
preimudata=LinkedList();% record the previous imu data

% read in imu data
[fimu, imudata, preimudata]=readimuheader(options.imufile, preimudata, options.startTime, Counter, imuFileType);
lastimu=preimudata.getLast();
preimutime=lastimu(1,end);
imuctr=1;   % to count how many IMU data after the latest GPS observations
% read the GPS data and align the GPS data with the imu data
if(useGPS)
    [fgps, gpsdata]=readgpsheader(gpsfile, preimutime, gpspostype);
else gpsdata=inf;
end
gpsctr=1;       % the number of GPS data that has been used

%load the camera frame timestamp data
if(options.useCam)
    % Set plot windows
    set_plots;
    % output video
    outputVideo = VideoWriter(videoName);
    outputVideo.FrameRate =20;
    open(outputVideo);
    camconfighistory=LinkedList();
    trajectory = zeros( 7, lastIm - initIm+1);
    generate_random_6D_sphere;
    [fimgtime, imgepoch, initIm]=readimgtimeheader(imgtimefile,preimutime, initIm);
    step=initIm;
    firstGroupFrameId=0;
    % Camera initialization
    [cam, resizeScale] = initialize_cam_v001(options.camtype, maxEdge);
    tracker= KLT_tracker(ftTblFile, sequencePath, imgseqtype, resizeScale);
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
    %Write IMU's position (actually antenna) and velocity, rpy and accel and gyro bias to the files
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
    if ((curimutime-covupt_time)>=options.maxCovStep) ||(curimutime>gpsdata(1)||(curimutime>imgepoch))
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
    isStatic =~isempty(options.zuptSE) && ~isempty(find(((options.zuptSE(:,1)<=curimutime)&(options.zuptSE(:,2)>=curimutime))==1,1));
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
            predict=filter.rvqs2e(4:6);
            H=sparse([zeros(3) skew(filter.rvqs2e(4:6)) zeros(3)...
                quat2dcm_v000(filter.rvqs2e(7:10)) zeros(3,size(filter.p_k_k,1)-12)]);
        end
        R=eye(3)*options.sigmaZUPT^2;
        filter.correctstates(predict,measure, H,R);
    end
    
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
    
    %Apply the  GPS observations
    if (curimutime>gpsdata(1))
        imuctr=0; % to count how many imu epochs after the recent gps observations
        measure=ecef2geo_v000([gpsdata(2)/180*pi;gpsdata(3)/180*pi;gpsdata(4)],1);        
        if(strcmp(filter.tag, 'EKF_IMU_GPS_EFRM'))
            lever=quatrot_v000(filter.rvqs2e(7:10),filter.Tant2imu,0);
            predict=filter.rvqs2e(1:3)+lever;
            H=[eye(3) zeros(3) skew(lever) zeros(3,size(filter.p_k_k,1)-9)];
        elseif(strcmp(filter.tag, 'EKF_IMU_GPS_NFRM'))
            gpsecef=ecef2geo_v000([gpsdata(2)/180*pi;gpsdata(3)/180*pi;gpsdata(4)],1);
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
        else
            R=4*diag(gpsdata(7:9).^2);
        end
        filter.correctstates(predict,measure, H,R);        
        %Read the next gps data
        [fgps, gpsdata]=grabnextgpsdata(fgps, gpspostype);
        gpsctr=gpsctr+1;
        if (gpsctr>=options.gpsnum)
            disp(['GPS outage starts from ' num2str(gpsdata(1)) 'GTOW sec!']);
            gpsdata(1)=inf;% stop using GPS, dead reckoning
        end
    end
    % image measurement
    if(curimutime>imgepoch)
        % there are some gap in the imu recordings
        if(abs(curimutime-imgepoch)>options.dt)
            disp(['abs(curimutime-imgepoch)>options.dt at ' num2str(curimutime)...
                ' GTOW sec and image index ' num2str(step) '!']);
        end
        % predict the image cooridnates of points in the states
        if(strcmp(filter.tag,'EKF_CAM_IMU_GPS'))
            filter.predictPoints(cam);
            [trackedPoints, pointGroup]=tracker.TrackandDetectFeaturePoints(...
                filter.features_info,filter.groupPose,filter.rvqs0([7:10, 1:3]),filter.camPose, cam, step, imgepoch );
        else
            [trackedPoints, pointGroup]=tracker.TrackandDetectFeaturePoints(...
                [],[],[filter.rvqs2e(7);-filter.rvqs2e(8:10); filter.rvqs2e(1:3)] ,filter.camPose, cam, step, imgepoch);
            assert(isempty(trackedPoints));
        end
        if(~isempty(trackedPoints))
            filter.SetFeaturesInfo(trackedPoints);
            % 1-Point RANSAC hypothesis and selection of low-innovation inliers
            filter.ransac_hypotheses( cam );
            % Partial update using low-innovation inliers
            filter.ekf_update_li_inliers( );
            filter.UpdateGrpPtsNo(step); % update the number of points in a group
            % "Rescue" high-innovation inliers
            filter.rescue_hi_inliers(  cam );
            % Partial update using high-innovation inliers
            %             filter.SetAllHiInliers();
            filter.ekf_update_hi_inliers();
            %             filter.ekf_update_rho( filter);
            %             disp([step, imgepoch, curimutime, length(filter.features_info)]);
            %             if(filter.camConfigSIP~=filter.groupFrameSIP)
            %                 if(camconfighistory.size()== Counter.numcamconfigrecords)
            %                     camconfighistory.removeFirst();
            %                 end
            %                 camconfighistory.addLast(filter.camPose);%record previous camera config estimate
            %                 % disable camera configuration, if they ramain almost constant over the
            %                 % past period, say 30 seconds
            %                 if(mod(step,5)==0&&(camconfighistory.size()==Counter.numcamconfigrecords))
            %                     [decision, estimate]=iscamerastable(listcontent(camconfighistory), 1*pi/180 ,0.03);
            %                     if(decision)
            %                         filter.disable_camerastates( estimate );
            %                         disp(['Camera boresight and bearing estimation disabled from ' num2str(curimutime) 'GTOW sec!']);
            %                     end
            %                 end
            %             end
            tracker.PostStates(filter.GetPtsDepth(), filter.rvqs0, step);
        end
          
        if(strcmp(filter.tag,'EKF_CAM_IMU_GPS'))
            plotFtPts;
            % (adding and deleting features
            filter.Renew( pointGroup);
            filter.SaveCamPoseandRqs02e(imgepoch,step);
            filter.LoopClosure();
            filter.CheckConstraints();
        elseif(~isempty(pointGroup))
            % recast the filter to EKF CAM IMU GPS
            if(strcmp(filter.tag,'EKF_IMU_GPS_NFRM'))             
                filter.SetCamAndRvqs2e(options);
            end
            filter=EKF_filter_s0frame(options, filter, pointGroup);
            % convert the e frame coordination of rs in e and qs2e, to s0
            % frame formulation, rs in s0, q s0 2s.
            tracker.FeatureTable_e2s0(filter.rqs02e);
            firstGroupFrameId=pointGroup(1).initFrmNo;
%             filter.disable_camerastates([]); % fix the camera states from the start
            filter.SaveCamPoseandRqs02e(imgepoch,step);
        end
        % Save images
        imgplot = getframe(gcf);
        writeVideo(outputVideo,imgplot);
        % take the next image time
        if(step<lastIm)
            hstream= fgetl(fimgtime);
            camel=sscanf(hstream,'%f,%f');
            imgepoch=camel(2);
            step=step+1;
            assert( camel(1)==step);
            % right now camel(1, end) is the time of this step image
        else
            imgepoch=inf; % stop using images
            close(outputVideo);
        end
    end
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
implay(videoName,30);
