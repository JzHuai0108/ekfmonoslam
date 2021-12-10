% This program tests ekf_filter_s0frame_bias class
% This script aims at low grade IMU dead reckoning 

% (1) microstrain IMU readings may have gaps as large as 20 secs, cut
% these, otherwise, gaps a bit larger than sampling interval are fine
% (2) a UAV has a collision upon landing, so ZUPT after landing is
% meaningless, so we discard data following and just prior to landing
% (3) Estimating states with ZUPT using static IMU data prior
% to launching seems problematic. I believe, it is due to accelerometer
% bias on vertical axis suddenly shoot up upon taking off. So it is advised
% to filter IMU and GPS data for UAV tests since launching
% (4) when a UAV flys up, it does not rotate itself until it reaches a
% certain height. Ideally in autonomous flight, a UAV only turns at waypoints
% (5) the UAV flies at maximum height of 25 meters, so we assume constant
% gravity
% (6) the UAV does not fly very far, so we assume fixed n-frame, and it
% does not fly very long, so earth rotation is ignored

% The experiment on 2014/10/08 collected data with Espon IMU with GPS,
% Microstrain 3DM-GX3-35 with GPS capability, 
% velodyne and Nikon cameras which were mounted on an Octorotor. But I
% have no access to these data at the moment, so the related tests are
% removed.

function uavmems_EKF_filter()
workspace_path = '/media/jhuai/docker/ekfmonoslam/';
addpath([workspace_path 'instk']); % imu functions
addpath([workspace_path 'voicebox']); % rotation functions
addpath([workspace_path 'ekfmonocularslamv02']); % filter classes
addpath([workspace_path 'utilities']); % data readers

% than dcm2quat_v000 and dcm2euler_v000
import java.util.LinkedList
clear variables;
clc; close all; format longg;
fprintf('\n EKF of MEMS IMU on UAV!\n\n');
rng('default');

experim=3;
switch experim
    case 1
        % test on Microstrain 3dm gx3-35 data integration with GPS data
        % collected by octoptor on 2014/10/08 velodyne test       
        isOutNED=true;
        resdir='C:\Users\huai.3\Desktop\huai work\OctoptorINSGPStest\temp\';
        filresfile=[resdir, 'filresult.bin']; % navigation states
        imuresfile=[resdir, 'imuresult.bin']; % imu error terms
        % imu options                       
        options.startTime=320334.7513; % for velodyne test just after launching
        options.endTime= 320522.6513; % exclude landing for velodyne test
        options.imuErrorModel=5; % 4 for random constant acc bias and gyro bias, 5 for random walk acc bias and random constant gyro bias
        options.mechanization=2; % 1 for wander azimuth
        
        options.imutype=5; % 5 for 3dm gx3 35, 7 for epson m-g362pdc1
        options.dt=1/100;  %sampling interval
        options.maxCovStep=options.dt; %maximum covariance propagation step, if equal to dt, means single speed mode
        options.imufile='F:\oct082014velodyne\IMU_MicroStrain-Velodyne.csv';
        
        imuFileType=4; % 4 for 3DM GX3 -35 csv, 5 for m-g362pdc1 csv
        inixyz_ant= [ 579144.6076            -4851437.6245           4086499.4477]';% the position of antenna at startTime for velodyne
        options.inillh_ant=ecef2geo_v000(inixyz_ant,0);
        Ce2n0=llh2dcm_v000(options.inillh_ant(1:2),[0;1]);
        
        options.imuErrors=zeros(6,1);         
        options.imuErrors(4:6)=[ 5.23585856565723e-07      0.000768720262789487     -3.76258753649856e-05]'; % for microstrain velodyne
        options.imuErrors(3)=1.5;
   
        options.Vn=[0;0;0]; % roughly estimated from GPS
        options.qb2n=  roteu2qr('xyz',[0;0; 80]/180*pi); % the rotation from imu to w-frame(rotated n-frame)
                
        options.InvalidateIMUerrors=false;
        options.initAttVar=2*pi/180; % 1 deg std for roll and pitch, 2 times 1 deg for yaw std
        
        % GPS options
        useGPS=false; % if to only use height observation, modify if(1) to if(0) in the following code
        useGPSstd=true; % use the std in the rtklib GPS solutons
        % for velodyne microstrain, useGPSstd gives better integration results
        Tant2imu=[8;0;-118]*1e-3; % from Greg's drawing, the position of antenna in microstrain frame
%         Tant2imu=[4;5;99]*1e-3; % from Greg's drawing, the position of antenna in epson frame
        % gps start and end time
        gpsSE=[options.startTime, options.startTime+2000];
        
        gpsfile='F:\oct082014velodyne\Velodyne.pos';
        isConstantVel=false; % use constant velocity model
        options.velNoiseStd=1; % velocity noise density m/s^2 in a horizontally axis
        
        % ZUPT start and end Time, n x 2 matrix, n is the number of segment
        zuptSE=[options.startTime, 320301.8429];%  320541.8429,  320688]; % for oct08 2014 velodyne dataset
        options.sigmaZUPT = 0.05;% unit m/s
        rateZUPT=round(sqrt(1/options.dt)); % from tests, we see zupt has adverse effect in this case
        % NHC options
        options.sigmaNHC = 0.1;% unit m/s
        rateNHC=inf; % turn off NHC
        % minimum velocity before applying the NHC, this option decouples ZUPT and NHC
        options.minNHCVel=2.0;
        
       % gravity correction related parameters
        rateGravity=inf;        
        gn0=-[ -0.180653756542627        0.0836986162201811     -9.80614572647246]'; % for microstrain velodyne test
        sigmaGravity=2; % unit m/s^2, gravity measurement noise
        gravityMagCutoff=2; 
        wie2n0=[0;0;0]; % no earth rotation
        % vertical velocity constraints
        rateVel=inf; %round(sqrt(1/options.dt)/2);
        sigmaVertVel=6; % unit m/s
        sigmaHorizVel=10;

    case 2
        % test on Epson data integration with GPS data
        % collected by GPS Van on 2015/11/11, session D
        isOutNED=true;
        resdir='C:\Users\huai.3\Desktop\huai work\OctoptorINSGPStest\temp\';
        filresfile=[resdir, 'filresult.bin']; % navigation states
        imuresfile=[resdir, 'imuresult.bin']; % imu error terms
        % imu options
        options.startTime=327674.0042;
        options.endTime= 327674.0042+ 1000; %329771.1481;

        options.imuErrorModel=5; % 4 for random constant acc bias and gyro bias, 5 for random walk acc bias and random constant gyro bias
        options.mechanization=2; % 1 for wander azimuth
        
        options.imutype=7; % 5 for 3dm gx3 35, 7 for epson m-g362pdc1
        options.dt=1/500;  %sampling interval
        options.maxCovStep=5*options.dt; %maximum covariance propagation step, if equal to dt, means single speed mode

        options.imufile='\\File.ceegs.ohio-state.edu\SPIN\MultiSensor_NOV_11_2015\IMU_Epson\20151111_140059_D_timetagged.csv';
        imuFileType=5; % 4 for 3DM GX3 -35 csv, 5 for m-g362pdc1 csv
        % the position of antenna at startTime
        inixyz_ant=[592574.6611  -4856604.0417   4078414.4645]';
        options.inillh_ant=ecef2geo_v000(inixyz_ant,0);
   
        options.imuErrors=zeros(6,1);        
        options.Vn=[0;0;0]; % roughly estimated from GPS
        options.Cb2imu=[1 0  0;  0 -1 0; 0 0 -1];% vehicle body frame to imu frame
        options.qb2n= roteu2qr('xyz', [0; 0; 8.9483]*pi/180); % obtained from H764G, yaw deg

        options.InvalidateIMUerrors=false;
        options.initAttVar=2*pi/180; % 1 deg std for roll and pitch, 2 times 1 deg for yaw std
        
        % GPS options
        useGPS=true; 
        useGPSstd=true; % use the std in the rtklib GPS solutons
        Tant2imu=[0.454; -0.746; 1.344];
        % gps start and end time
        gpsSE=options.startTime+[0, 350; 420, 500; 540, 600];
        
        gpsfile='C:\JianzhuHuai\GPS_IMU\programs\matlab_ws\data\Novatel_rover_KIN_20151111.pos';
        
        isConstantVel=false; % use constant velocity model
        options.velNoiseStd=1; % velocity noise density m/s^2 in a horizontally axis
        
        % ZUPT options
        % ZUPT start and end Time, n x 2 matrix, n is the number of segment        
        zuptSE=[options.startTime, options.startTime+180];
        options.sigmaZUPT = 0.1;% unit m/s
        rateZUPT=round(sqrt(1/options.dt));
        % NHC options
        options.sigmaNHC = 0.1;% unit m/s
        rateNHC=inf; % turn off NHC
        % minimum velocity before applying the NHC, this option decouples ZUPT and NHC
        options.minNHCVel=2.0;
        
        % gravity correction related parameters
        rateGravity=inf;
        gn0= [ 0; 0; 9.804]; % gravity in n frame at start
        sigmaGravity=2; % unit m/s^2, gravity measurement noise
        gravityMagCutoff=2; 
        
        wie2n0=[0;0;0]; % earth rotation
        
        % vertical velocity constraints
        rateVel=inf; %round(sqrt(1/options.dt)/2);
        sigmaVertVel=6; % unit m/s
        sigmaHorizVel=10; 
    case 3
        % test on inertial data captured by the Epson IMU mounted on a car.
        % The GPS is kept out in the last 20 seconds to verify that the
        % INS works in a short span.
        isOutNED=true;
        resdir='/home/jhuai/Desktop/temp/gnssimu/';
        filresfile=[resdir, 'filresult.bin']; % navigation states
        imuresfile=[resdir, 'imuresult.bin']; % imu error terms
        % imu options
        options.startTime=327674.0042;
        options.endTime= 327674.0042+ 280;

        options.imuErrorModel=5; % 4 for random constant acc bias and gyro bias, 5 for random walk acc bias and random constant gyro bias
        options.mechanization=2; % 1 for wander azimuth
        
        options.imutype=7; % 5 for 3dm gx3 35, 7 for epson m-g362pdc1
        options.dt=1/30;  %sampling interval
        options.maxCovStep=options.dt; %maximum covariance propagation step, if equal to dt, means single speed mode

        options.imufile='/media/jhuai/SeagateData/jhuai/data/osu-spin-lab/20151111/20151111_140059_D_timetagged.csv';
        imuFileType=5; % 4 for 3DM GX3 -35 csv, 5 for m-g362pdc1 csv
        % the position of antenna at startTime
        inixyz_ant=[592574.6611  -4856604.0417   4078414.4645]';
        options.inillh_ant=ecef2geo_v000(inixyz_ant,0);
   
        options.imuErrors=zeros(6,1);        
        options.Vn=[0;0;0]; % roughly estimated from GPS
        options.Cb2imu=[1 0  0;  0 1 0; 0 0 1];% vehicle body frame to imu frame
        accel= [0 0 9.81]; % The average accelerometer measurement vector at the beginning.
        eulers2n=orientbymagnaccel(accel);
        options.qb2n= roteu2qr('xyz', eulers2n); % estimated from accelerometer

        options.InvalidateIMUerrors=false;
        options.initAttVar=2*pi/180; % 1 deg std for roll and pitch, 2 times 1 deg for yaw std
        
        % GPS options
        useGPS=true; 
        useGPSstd=true; % use the std in the rtklib GPS solutons
        Tant2imu=[0.454; -0.746; 1.344];
        % gps start and end time
        gpsSE=options.startTime+[0, 260];
        
        gpsfile='/media/jhuai/SeagateData/jhuai/data/osu-spin-lab/20151111/Novatel_rover_KIN_20151111.pos';
        
        isConstantVel=false; % use constant velocity model
        options.velNoiseStd=1; % velocity noise density m/s^2 in a horizontally axis
        
        % ZUPT options
        % ZUPT start and end Time, n x 2 matrix, n is the number of segment        
        zuptSE=[options.startTime, options.startTime+30];
        options.sigmaZUPT = 0.1;% unit m/s
        rateZUPT=round(sqrt(1/options.dt));
        % NHC options
        options.sigmaNHC = 0.1;% unit m/s
        rateNHC=inf; % turn off NHC
        % minimum velocity before applying the NHC, this option decouples ZUPT and NHC
        options.minNHCVel=2.0;
        
        % gravity correction related parameters
        rateGravity=inf;
        gn0= [ 0; 0; 9.804]; % gravity in n frame at start
        sigmaGravity=2; % unit m/s^2, gravity measurement noise
        gravityMagCutoff=2; 
        
        wie2n0=[0;0;0]; % earth rotation
        
        % vertical velocity constraints
        rateVel=inf; %round(sqrt(1/options.dt)/2);
        sigmaVertVel=6; % unit m/s
        sigmaHorizVel=10; 
    otherwise
        error(['Unsupported testing case ' num2str(experim) '!']);
end
numPrevImuDataToKeep=40; % the number of data put in preimudata
% Initialize the model state and covariance of state, process noise and
% measurment noise

filter =EKF_filter_s0frame_bias(options);

preimudata=LinkedList();% record the previous imu data

% read in imu data
[fimu, imudata, preimudata]=readimuheader(options.imufile, preimudata, options.startTime, numPrevImuDataToKeep, imuFileType);
lastimu=preimudata.getLast();
preimutime=lastimu(1,end);

gpsCount=0;
imuCountSinceGnss=1;   % to count how many IMU data after the latest GPS observations
% read the GPS data and align the GPS data with the imu data
if(useGPS)
    [fgps, gpsdata, gpspostype]=readgpsheader(gpsfile, gpsSE(1,1), 'lla');
else
    gpsdata=inf;
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
    %% Write IMU's position and velocity, rpy and accel and gyro bias to the files
    if(isOutNED)
        filter.SaveToFile(options.inillh_ant, preimutime, ffilres);
    else
        filter.SaveToFile([], preimutime, ffilres);
    end
    if (imudata(1)-initime)>60
        disp(['Process Time:' num2str(imudata(1))]);
        initime=imudata(1);
    end
    % state propagation
    imuaccum=filter.ffun_state(imuaccum,[imudata(2:7);preimutime;curimutime; gn0; wie2n0], 0, isConstantVel);
    
    % covariance propagation
    if ((curimutime-covupt_time)>=options.maxCovStep)
        %propagate the covariance
        filter.ffun_covariance(imuaccum, covupt_time, curimutime, isConstantVel);
        covupt_time=curimutime;
        imuaccum=zeros(6,1);
        %Record covariance and navigation solution for the smoother
        %Note that I am recording the predicted solutions. That is why the
        %backward part must use filtered solution.
        
        fwrite(fimures,[curimutime;filter.imuErrors;...
            full(sqrt(diag(filter.p_k_k(filter.imuBiasDriftSIP+(0:5),...
            filter.imuBiasDriftSIP+(0:5)))))],'double');
    end
    %% ZUPT (zero velocity updates).
    isStatic =~isempty(zuptSE) && ~isempty(find(((zuptSE(:,1)<=curimutime)&(zuptSE(:,2)>=curimutime))==1,1));
    isZUPT =mod(imuCountSinceGnss, rateZUPT)==0;
    if (isStatic&&isZUPT)
        measure=zeros(3,1);
        predict=filter.rvqs0(4:6);
        H=sparse([zeros(3) eye(3) zeros(3,size(filter.p_k_k,1)-6)]);
        R=eye(3)*options.sigmaZUPT^2;
        filter.correctstates(predict,measure, H,R);
    end
    %% Gravity norm measurements.
    useGravity =mod(imuCountSinceGnss, rateGravity)==0;
    if (useGravity&& abs(norm(imudata(2:4))-9.81)< gravityMagCutoff)
        predAccel=imudata(2:4)-filter.imuErrors(1:3);
        predict=-quatrot_v000(filter.rvqs0(7:10),predAccel,1);
        H=[zeros(3,6) rotqr2ro(filter.rvqs0(7:10))'*skew(predAccel) zeros(3,filter.covDim-9)];
        measure=gn0;
        R=sigmaGravity^2*eye(3);
        filter.correctstates(predict,measure, H,R);
    end
    %% NonHolonomic constraints.
    isNHC =mod(imuCountSinceGnss, rateNHC)==0;
    velnorm=filter.GetVelocityMag();
    if (isNHC&&velnorm>options.minNHCVel)
        % non-holonomic constraints
        Cs02s=quat2dcm_v000(filter.rvqs0(7:10));
        Cs02b=filter.Cb2imu'*Cs02s;
        predict=Cs02b(2:3,:)*filter.rvqs0(4:6);
        larry=filter.Cb2imu'*skew(quatrot_v000(filter.rvqs0(7:10),filter.rvqs0(4:6),0));
        H=sparse([zeros(2,9) Cs02b(2:3,:) larry(2:3,:) zeros(2,size(filter.p_k_k,1)-15)]);
        
        measure=[0;0];
        R=eye(2)*options.sigmaNHC^2;
        filter.correctstates(predict,measure, H,R);
    end
    %% Vertical velocity constraints.
    useVel =mod(imuCountSinceGnss, rateVel)==0;
    if (useVel)        
        predict=filter.rvqs0(4:6);
        H=[eye(3), zeros(3,filter.covDim-3)];
        measure=zeros(3,1);
        R=diag([sigmaHorizVel,sigmaHorizVel,sigmaVertVel].^2);
        filter.correctstates(predict,measure, H,R, 1);
    end
    
    %% GNSS observations.
    if (curimutime>=gpsdata(1))
        gpsCount = gpsCount + 1; 
        if rem(gpsCount, 20) == 0 
            fprintf('Using GNSS data at %.3f.\n', gpsdata(1, 1)); 
        end
        imuCountSinceGnss=0;
        % get antenna position in N frame.
        % 1. matlab mapping toolbox approach
        [east, north, up] = geodetic2enu(gpsdata(2), gpsdata(3), gpsdata(4), ...
            options.inillh_ant(1), options.inillh_ant(2), options.inillh_ant(3), ...
            wgs84Ellipsoid, 'radians');
        p_N_ant = [north; east; -up];
        % 2. our own implementation
%         gpsecef=lla2ecef([gpsdata(2:3) * 180 / pi; gpsdata(4)]')';
%         p_N_ant_goodold = quatrot_v000(filter.rqs02e(4:7), gpsecef - inixyz_ant, 1);
%         assert(max(abs(p_N_ant - p_N_ant_goodold)) < 1e-8);
        % 3. matlab automated driving toolbox approach results differs
        % from the first two approaches by 2 mm.
%         [east2, north2, up2] = latlon2local(gpsdata(2), gpsdata(3), gpsdata(4), options.inillh_ant);

        measure = p_N_ant - quatrot_v000(filter.rvqs0(7:10),Tant2imu,1);

        predict=filter.rvqs0(1:3);
        % use three channels
        if(1)
            H=sparse([eye(3), zeros(3), zeros(3,size(filter.p_k_k,1)-6)]);
            % the following setting of noise variances is suitable for RTKlib output
            if(~useGPSstd)
                if(gpsdata(5)==1)
                    R=diag([0.05,0.05,0.1].^2);
                elseif(gpsdata(5)==2)
                    R=diag([0.5,0.5,1.0].^2);
                else
                    R=diag([10,10,10].^2);
                end
            else
                R=4*diag(gpsdata(7:9).^2);
            end
            filter.correctstates(predict,measure, H,R);
        else % only use height
            H=sparse([0,0,1, zeros(1,size(filter.p_k_k,1)-3)]);
            % the following setting of noise variances is suitable for RTKlib output
            if(~useGPSstd)
                if(gpsdata(5)==1)
                    R=0.05^2;
                elseif(gpsdata(5)==2)
                    R=0.1^2;
                else
                    R=1^2;
                end
            else
                R=4*gpsdata(9)^2;
            end
            filter.correctstates(predict(3),measure(3), H,R, 1);
        end
        %Read the next gps data that is within the specified sessions
        gpsSErow=find(((gpsSE(:,1)<=gpsdata(1))&(gpsSE(:,2)>=gpsdata(1)))==1,1); % on which row/ session is the last gpsdata
        lastgpstime =gpsdata(1);             
        assert(~isempty(gpsSErow));
        [fgps, gpsdata]=grabnextgpsdata(fgps, gpspostype, 'lla');
        if (gpsdata(1)>gpsSE(gpsSErow,2))
            disp(['GPS outage starts from ' num2str(lastgpstime) ...
                ' GTOW sec which is ', num2str(lastgpstime - options.startTime), ...
                ' since the start!']);
            % find the next session
            nextSErow=0;
            gpsSErow=find(((gpsSE(:,1)<=gpsdata(1))&(gpsSE(:,2)>=gpsdata(1)))==1,1);
            if(isempty(gpsSErow))
                for iota= 1: size(gpsSE,1)
                    if gpsSE(iota,1)>= gpsdata(1)
                        nextSErow= iota;                        
                        break;
                    end
                end
                if(nextSErow==0)
                    disp(['GPS completely gone from ' num2str(lastgpstime) ...
                          ' GTOW sec which is ', num2str(lastgpstime - options.startTime), ...
                          ' since the start!']);
                    gpsdata(1)=inf;% completely stop using GPS
                else
                    while(gpsdata(1)< gpsSE(nextSErow,1))
                        [fgps, gpsdata]=grabnextgpsdata(fgps, gpspostype);                        
                        while(gpsdata(1) > gpsSE(nextSErow,2) && nextSErow< size(gpsSE,1))
                            nextSErow= nextSErow+1;                            
                        end
                    end
                    disp(['GPS will resume from ' num2str(gpsdata(1)) ...
                          ' GTOW sec which is ', num2str(gpsdata(1) - options.startTime), ...
                          ' since the start!']);            
                end
            else
                nextSErow = gpsSErow;
                disp(['GPS resumes from ' num2str(gpsdata(1)) ...
                      ' GTOW sec which is ', num2str(gpsdata(1) - options.startTime), ...
                      ' since the start!']);
            end         
        end
    end
    
    % Read the next imu data
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
        imuCountSinceGnss=imuCountSinceGnss+1;
    end
end
fprintf('Done.\n\n');
fclose all;

%% Plot navigation results
kf = readdata(filresfile, 1+18);

% Extracting GPS data
posdata=loadAllGPSData(gpsfile, [options.startTime, options.endTime]); % reads specified columns
% in given example these columns refer to: TOW, X, Y, Z, sol type, mX, mY, Mz
inillh_ant=ecef2geo_v000(inixyz_ant,0);
Ce2n0=llh2dcm_v000(inillh_ant(1:2),[0;1]);
qs02e=rotro2qr(Ce2n0'); % assume body frame is sensor frame
xs0=zeros(size(posdata,1),3);
for i=1:size(posdata,1)
    xs0(i,:)=quatrot_v000(qs02e, posdata(i,2:4)'-inixyz_ant ,1)';
end

nextFig=1;
f(nextFig) = figure;
plot3(kf(:,1)-kf(1,1),kf(:,3),kf(:,2),'g.');
hold on;
plot3(posdata(:,1)-kf(1,1),xs0(:,2),xs0(:,1),'r+');
grid;
axis equal;
xlabel('Time [s]');
ylabel('East [m]');
zlabel('North[m]');
legend('Trajectory of the IMU sensor per EKF', 'Trajectory of the antenna per GNSS');
title('Planar trajectory over time');

nextFig=nextFig+1;
f(nextFig) = figure;
plot3(kf(:,2),kf(:,3), kf(:, 4), 'g.');
hold on;
plot3(xs0(:,1),xs0(:,2), xs0(:, 3), 'r+');
grid;
axis equal;
xlabel('North [m]');
ylabel('East [m]');
zlabel('Down [m]');
legend('Trajectory of the IMU sensor per EKF', 'Trajectory of the antenna per GNSS');
title('Trajectory');

nextFig=nextFig+1;
f(nextFig) = figure;
plot(kf(:,1)-kf(1,1),kf(:,4),'g.')
hold on
plot(posdata(:,1)-kf(1,1),xs0(:,3),'r+');

grid
xlabel('Time [s]')
ylabel('Height/m')
title('Height of antenna by KF(green) and reference(red)');
saveas(f(nextFig),[resdir 'red truth and height'],'fig');

filresfile=[resdir, 'filresult.bin'];
kf = readdata(filresfile, 1+18);
plotkf_v001(kf, resdir);
imuresfile=[resdir, 'imuresult.bin'];
err = readdata(imuresfile, 1+12);

nextFig=nextFig+1;
nextDim=2;
f(nextFig) = figure;
loc = nextDim+(0:2);
plot(err(:,1)-err(1,1),err(:,loc),'marker','.')
grid
xlabel('Time [s]')
ylabel('m/s^2')
legend('X','Y','Z')
title('Accelerometer bias drift');

%close all
nextFig=nextFig+1;
nextDim=nextDim+3;
f(nextFig) = figure;
loc = nextDim+(0:2);
plot(err(:,1)-err(1,1),err(:,loc),'marker','.')
grid
xlabel('Time [s]')
ylabel('radian/s')
legend('X','Y','Z')
title('Gyro bias drift');

%close all
nextFig=nextFig+1;
nextDim=nextDim+3;
f(nextFig) = figure;
loc = nextDim+(0:2);
plot(err(:,1)-err(1,1),err(:,loc),'marker','.')
grid
xlabel('Time [s]')
legend('X','Y','Z')
title('Accelerometer bias drift cov std');

%close all
nextFig=nextFig+1;
nextDim=nextDim+3;
f(nextFig) = figure;
loc = nextDim+(0:2);
plot(err(:,1)-err(1,1),err(:,loc),'marker','.')
grid
xlabel('Time [s]')
legend('X','Y','Z')
title('Gyro bias drift cov std');

end

