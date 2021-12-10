% This program implements the Kinect IMU integration as discussed by Huai
% Formulated in a local frame with camera measurements
% the output trajectory is that of the IMU

% we face the problem that some IMU records are missing, but it is
% generally very small gap, so it is reasonable to use the received
% measurements for Kalman update once the imu epoch exceeds their epoch.

% Test cases (1) start IMU and Kinect in stationary mode, do ZUPT for IMU
% for at least 1 minutes, the initial frame of IMU and the camera are s0
% and c0 frame, then when Kinect is moved smoothly in predefined
% trajectory, the filter is used to estimate gravity in s0 frame, the
% calibration terms of camera and IMU

function imu_cam_calib()
addpath('..\instk'); % imu functions
addpath('..\voicebox\'); % for rotro2qr and rotqr2eu, they are more robuts
% than dcm2quat_v000 and dcm2euler_v000
import java.util.LinkedList
clear variables;
clc; close all; format longg;
fprintf('\n IMU_Kinect_Calib.m to test EKF filtering in calibrating mems IMU and Kinect!\n\n');
rng('default');

experim=4; %test case for IMU CAMera calibration

switch experim
    
    case 4
        % test on Steval MKI062V2 data integration with camera, Mar 2014
        % output options
        % test observations: block camera measurement in the ZUPT period,
        % or turn off ZUPT, does not change the results much of turn-on ZUPT and camera
        % measurement at the same time.
        % the strange profile of Cs2c has to do with the data,
        % Ts2c wrong becuase of the filter mechanism and data.
        % use of gravity magnitude constraint does not show much improvement
        % it is advised to turn off scale factor estimate, let gravity and
        % RS2C float
        isOutNED=true;
        resdir='C:\Users\huai.3\kinectfusion\KinectFusionExplorer-D2D\temp\';
        filresfile=[resdir, 'filresult.bin']; % navigation states
        imuresfile=[resdir, 'imuresult.bin']; % imu error terms        
        % imu options
        options.startTime=2;
        options.endTime=130;
        options.imuErrorModel=4; % how bias and scale factor is modeled
        options.mechanization=2; % 1 for wander azimuth       
        options.imutype=6;      %  Steval MKI062V2
        options.dt=1/50;  %sampling interval
        options.maxCovStep=options.dt; %maximum covariance propagation step, if equal to dt, means single speed mode
        options.Cb2imu=eye(3); 
        options.Timu2body=zeros(3,1); % h764G is the body frame
        options.imufile=[resdir,'iNemo_02092014_061806.tsv'];
        
       % Toc2s=[[  0.0970;   -0.2459;    0.9432];-[   -0.0037;   -0.0537;   -0.0680]]; % PoseIMUKinectx302010
       % Toc2s=[[   0.0886;   -0.2499;    0.9400]; -[  -0.0006;   -0.0526;   -0.0681]]; % for PoseIMUKinecty102030
      %  Toc2s=[[0.1194;   -0.2479;    0.9359]; -[  -0.0022;   -0.0530;   -0.0676]];% for PoseIMUKinectz201306
      %  Toc2s=[[   0.177695680320824;-0.00680201661415064;0.962107115439917]; -[  0.0362812169298294;-0.0710746054253426;0.0478464512175324]]; % for PoseIMUCircle
       % Toc2s=[[   -0.219889874055416; -0.230975516372796;0.907607304785909]; -[ -0.0084376;   -0.057877;       -0.068204]]; % for PoseIMUKinect3
       Toc2s=[0.285567755532141, -0.0225589251844046, 1.0874587565859,  0, 3/180*pi, 4/180*pi]'; % for iNEMO2_20140603_160228 anticlock
        % ZUPT start and end Time, determine before the filter
        % format n x 2
        % n is the number of segment
        % first  column is the start time
        % second column is the end time, for example [575908      576243; 576502      576602]        
        options.zuptSE=[2, 50]; % for iNEMO2_20140603_160228 anticlock
        imuFileType=3; % iNemo data
       % imuFileType=2; % file provided by Yujia, mix of steval and kinect output
        %Initial PVA of IMU
        options.inillh_ant=[40.00311241687*pi/180; -83.01529294250*pi/180; 227.901];
        % lat, log, height of the antenna, in this case, identical to IMU,
        % guessed from http://www.daftlogic.com/sandbox-google-maps-find-altitude.htm
        options.imuErrors=zeros(12,1);       
        options.imuErrors(1:6)=Toc2s; % for the Steval case, comment this line if not steval
        options.Vn=[0;0;0];
        options.Ve=[0;0;0];
        options.qb2n=rotro2qr(R3(0.1477)*[0, -1,0; -1,0, 0; 0,0, -1]);
        % 0.1477 is roughly estimated from earth surface triangulation using data from daftlogic
        options.InvalidateIMUerrors=false;
        options.initAttVar=0.5*pi/180; % 1 deg std for roll and pitch, 2 times 5 deg for yaw std
        % gps options
        useGPS=false;
        
        % ZUPT options
        options.sigmaZUPTPos = 0.02;% unit m
        options.sigmaZUPTVel = 0.02;% unit m/s
        options.sigmaZUPTAng = 1*180/pi;% unit rad
        rateZUPT=round(sqrt(1/options.dt)); % from tests, we see zupt has adverse effect in this case
        % NHC options
        options.sigmaNHC = 0.1;% unit m/s
        rateNHC=inf; %round(sqrt(1/options.dt));
        % minimum velocity before applying the NHC, this option decouples ZUPT and NHC
        options.minNHCVel=2.0;
        
        rateGravNorm=inf; %rateZUPT;
        sigmaGravMag=3e-3; % unit m/s^2
        
        options.useCam=false;
        camFile=options.imufile; % where the camera measurements come from
        options.camPoseFile=[resdir, 'kinectPose.txt']; % output camera position and attitude
        deltaPhi=[0;0;0]; % [-0.487; 1.106; -0.518]/180*pi;
        options.Cimu2cam=rotqr2ro(rvec2quat_v000(deltaPhi))*[0,0,1;-1,0,0;0,-1,0]'; % (Rc2s)'
        options.Tcam2body=[3; -6.2; -4.05]*1e-2; % unit m
        % depth camera in body frame whichStep is assumed to be the imu frame,
        
        sigmaEuler=[1; 1; 1]*pi/180; % unit rad
        sigmaTrans=[2;2;2]*1e-3; % unit m
    otherwise
        error('Unsupported testing case!');
end
numPrevImuDataToKeep=40; % the number of data put in preimudata
% Initialize the model state and covariance of state, process noise and
% measurment noise
filter =EKF_filter_s0framelet(options);
preimudata=LinkedList();% record the previous imu data

% read in imu data
[fimu, imudata, preimudata]=readimuheader(options.imufile, preimudata, options.startTime, numPrevImuDataToKeep, imuFileType);
lastimu=preimudata.getLast();
preimutime=lastimu(1,end);
imuctr=1;   % to count how many IMU data after the latest GPS observations

%load the camera measurements with timestamp
if(options.useCam)
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
    filter.SaveToFile([], preimutime, ffilres);
    if (imudata(1)-initime)>60
        disp(['Process Time:' num2str(imudata(1))]);
        initime=imudata(1);
    end
    % time propagation of the IMU
    imuaccum=filter.ffun_state(imuaccum,[imudata(2:7);preimutime;curimutime]);
    
    %Update the covariance
    if (((curimutime-covupt_time)>=options.maxCovStep) ||(curimutime>=imgepoch))
        %propagate the covariance
        filter.ffun_covariance(imuaccum, covupt_time, curimutime);
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
        measure=zeros(9,1);
        predict=[filter.rvqs0(1:6); unskew(eye(3)-rotqr2ro(filter.rvqs0(7:10)))];
        H=sparse([eye(9) zeros(9,size(filter.p_k_k,1)-9)]);
        R=zeros(9); 
        R(1:3,1:3)=eye(3)*options.sigmaZUPTPos^2;
        R(4:6,4:6)=eye(3)*options.sigmaZUPTVel^2;
        R(7:9,7:9)=eye(3)*options.sigmaZUPTAng^2;
        filter.correctstates(predict,measure, H,R);
    end
    % apply gravity norm measurement
    isGravNorm =mod(imuctr, rateGravNorm)==0;
    if (isStatic&&isGravNorm)
        %Gravity (most time consuming part of ecef implementations)
        xyz_imu=filter.rqs02e(1:3)+quatrot_v000(filter.rqs02e(4:7),filter.rvqs0(1:3),0);
        Llh=ecef2geo_v000(xyz_imu,0);
        [Rn, Re, gn, sL, cL, WIE_E]=geoparam_v000(Llh);
        predict=sqrt(imudata(2:4)'*imudata(2:4));
        gs0=quatrot_v000(filter.rvqs0(7:10), imudata(2:4),1); % predicted gravity in s0 frame
        H=[zeros(1,9) 1/predict*gs0' zeros(1,filter.covDim-12)];
        measure=norm(gn,2);
        R=sigmaGravMag^2;
        filter.correctstates(predict,measure, H,R);
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
    if(preimudata.size()==numPrevImuDataToKeep)
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
% draw Rc2c0 predicted and measured values
% nextFig=1;
% f(nextFig) = figure;
% plot(camdata(:,1)-camdata(1,1),camdata(:, 2:4)*180/pi,'marker','+');
% hold on
% plot(euc2c0imu(:,1)-euc2c0imu(1,1),euc2c0imu(:, 2:4)*180/pi,'marker','v');
% grid
% xlabel('Time [s]')
% ylabel('degree')
% legend('Roll cam','Pitch cam','Yaw/Heading cam','Roll imu','Pitch imu','Yaw/Heading imu');
% title('Cc2c0 by Kinect and predicted by IMU');

fclose all;
% display the navigation results
plotnav;
end

function RoTrKinect=load_KinectDataset(dataFile, startTime)
% output time, euler angle of Rc2c0, T c0 in c
fid = load(dataFile); %'PoseIMUKinect3.tsv');
DataIMU = fid((1:2:size(fid,1)-1),:);
EulerKin = fid((2:2:size(fid,1)),2:4);
TranKin = fid((2:2:size(fid,1)),5:7);
obsindex = zeros(size(EulerKin,1),1);
for i = 2:size(EulerKin,1)
    if EulerKin(i,1) ~= EulerKin(i-1,1)
        obsindex(i) = 1;
    end
end
whom=logical(obsindex)& (DataIMU(:,1)>=startTime);
RoTrKinect=[DataIMU(whom,1), EulerKin(whom, :), TranKin(whom, :)];
end
