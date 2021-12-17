% This program tests ekf_filter_s0frame_bias class which is suitable for
% GNSS + low-grade IMU fusion.

% Test findings:
% (1) Microstrain 3DM-GX3-35 IMU readings may have gaps as large as 20 secs.
% (2) The UAV has a collision upon landing which changes the biases
% dramatically.
% (3) The estimated biases when the car is static but with engine running, 
% or when the UAV is about to launch may be not so reliable as default zeros. 
% (4) when a UAV flys up, it does not rotate itself until it reaches a
% certain height. Ideally in autonomous flight, a UAV only turns at waypoints
% (5) the UAV flies at a max height of dozens of meters, say 25 m, so we may assume constant
% local gravity.
% (6) the UAV does not fly very far, so we may ignore the earth curvature.
% For MEMS IMU, we can safely ignore the Coriolis forces.

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

clear variables;
clc; close all; format longg;
rng('default');

experim=3;
switch experim
    case 1
        % test with multiple GNSS and IMU data collected by sensors on a 
        % scout mini on 2021/11/27.
        isOutNED=true;
        resdir='/home/jhuai/Desktop/temp/gnssimu/';
        filresfile=[resdir, 'filresult.csv']; % navigation states
        imuresfile=[resdir, 'imuresult.csv']; % imu error terms
        % imu options
        options.imuErrorModel=5; % 4 for random constant acc bias and gyro bias, 5 for random walk acc bias and random constant gyro bias
        options.mechanization=2;
        
        options.imutype=6; % 6 for steval mki062v2
        options.dt=1/230;  %sampling interval
        options.maxCovStep= 4 * options.dt; %maximum covariance propagation step, if equal to dt, means single speed mode
        imufile='/media/jhuai/Seagate2TB/jhuai/data/gnss-imu/20211127RTK_IMU/imu/com34.txt'; 
        imutimefile='/media/jhuai/Seagate2TB/jhuai/data/gnss-imu/20211127RTK_IMU/imu/com34-syncedlocaltimes.txt'; 
        allImuData = readmatrix(imufile);
        allImuTimes = readmatrix(imutimefile);
        allImuData(:, 1) = allImuTimes;
        allGyroData = allImuData(:, 2:4);
        allImuData(:, 2:4) = allImuData(:, 5:7);
        allImuData(:, 5:7) = allGyroData;
        dt_imu_to_gps = 0.0;
        allImuData(:, 1) = allImuData(:, 1) + dt_imu_to_gps;

        options.inillh_ant = [30.5298813668333 * pi / 180, 114.350355321667 * pi / 180, 20.8869]';
        options.Cb2imu=[1 0  0; 0 1 0; 0 0 1];% vehicle body frame to imu frame 
        options.Vn=[0;0;0];
        R_ned_enu = [0, 1, 0; 1, 0, 0; 0, 0, -1]; 
        R_W_b = R_ned_enu * R3(-45 * pi / 180); 
        options.qb2n = rotro2qr(R_W_b);

        options.InvalidateIMUerrors=false;
        options.initAttVar=2*pi/180; % std for roll and pitch, 3 times std for yaw

        % GPS options
        useGPS=true;
        useGPSstd=true; % use the std provided by the GPS data.
        options.imu_p_ant = [3.5; 3.5; 11.2] * 1e-2; % the position of antenna in IMU frame

        isConstantVel=false; % use constant velocity model
        options.velNoiseStd=1; % velocity noise density m/s^2 in a horizontally axis
 
        gpsfile='/media/jhuai/Seagate2TB/jhuai/data/gnss-imu/20211127RTK_IMU/processed/ublox_20211127_lla.txt'; 
        allGpsData = readmatrix(gpsfile); 
        allGpsData(:, 2:3) = allGpsData(:, 2:3) * pi / 180; 
        allGpsData(:, 7:9) = allGpsData(:, 5:7); 
        allGpsData(:, 5:6) = 0;

        % align GPS and IMU data at the beginning 
        imubeginbuffer = options.dt * 5;
        maxstarttime = max(allGpsData(1, 1), allImuData(1, 1)) + imubeginbuffer;  
        allGpsData = allGpsData(allGpsData(:, 1) >= maxstarttime, :);
        allImuData = allImuData(allImuData(:, 1) >= maxstarttime - imubeginbuffer, :);

        options.startTime = maxstarttime;
        options.endTime = maxstarttime + 600;
        gpsSE=[options.startTime, options.endTime];
        gpsEndIndex = find(allGpsData(:, 1) > gpsSE(2), 1, 'first');
        allGpsData = allGpsData(1:gpsEndIndex, :);

        % ZUPT start and end Time, n x 2 matrix, n is the number of segment
        zuptSE=[options.startTime, options.startTime + 300];
        options.sigmaZUPT = 0.1; % unit m/s
        numImuDataPerZUPT=round(sqrt(1/options.dt));

        % NHC options
        options.sigmaNHC = 0.1;% unit m/s
        numImuDataPerNHC = round(sqrt(1/options.dt));
        % minimum velocity before applying the NHC, this option decouples ZUPT and NHC
        options.minNHCVel=1.0;

        gn0=[-0.00005; 0.00036; 9.79347]; % WHU xinghu building rear door
        wie2n0=[0;0;0]; % no earth rotation

        options.imuErrors=zeros(6,1);
        staticUntilIndex = find(allImuData(:, 1) < zuptSE(2), 1, 'last' );
        gyroBias = mean(allImuData(1:staticUntilIndex, 5:7), 1)';  
        accelBias = mean(allImuData(1:staticUntilIndex, 2:4), 1)' + ...
            options.Cb2imu * quatrot_v000(options.qb2n, gn0, 1);  
        options.imuErrors = [accelBias; gyroBias];
    case 2
        % test on Epson data integration with GPS data
        % collected by GPS Van on 2015/11/11, session D.
        isOutNED=true;
        resdir='/home/jhuai/Desktop/temp/gnssimu/';
        filresfile=[resdir, 'filresult.csv']; % navigation states
        imuresfile=[resdir, 'imuresult.csv']; % imu error terms
        % imu options
        options.startTime=327674.0042;
        options.endTime= 327674.0042 + 290;

        options.imuErrorModel=5; % 4 for random constant acc bias and gyro bias, 5 for random walk acc bias and random constant gyro bias
        options.mechanization=2;
        
        options.imutype=7; % 7 for epson m-g362pdc1
        options.dt=1/500;
        options.maxCovStep=5*options.dt; %maximum covariance propagation step, if equal to dt, means single speed mode
        imufile='/media/jhuai/SeagateData/jhuai/data/osu-spin-lab/20151111/20151111_140059_D_timetagged.csv';
        imuFileType=5; % 5 for m-g362pdc1 csv
        allImuData=loadImuData(imufile, imuFileType, [options.startTime - 0.1, options.endTime + options.dt]); 
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
        options.imu_p_ant = [0.454; -0.746; 1.344];
        % gps start and end time
        gpsSE=options.startTime+[0, 270];

        gpsfile='/media/jhuai/SeagateData/jhuai/data/osu-spin-lab/20151111/Novatel_rover_KIN_20151111.pos'; 
        allGpsData=loadAllGPSData(gpsfile, [options.startTime, options.endTime], 'lla'); 

        isConstantVel=false; % use constant velocity model
        options.velNoiseStd=1; % velocity noise density m/s^2 in a horizontally axis
        
        % ZUPT options
        % ZUPT start and end Time, n x 2 matrix, n is the number of segment        
        zuptSE=[options.startTime, options.startTime+180];
        options.sigmaZUPT = 0.1;% unit m/s
        numImuDataPerZUPT=round(sqrt(1/options.dt));
        % NHC options
        options.sigmaNHC = 0.1;% unit m/s
        numImuDataPerNHC=inf; % turn off NHC
        % minimum velocity before applying the NHC, this option decouples ZUPT and NHC
        options.minNHCVel=2.0;
        
        gn0= [ 0; 0; 9.804]; % gravity in n frame at start
        
        wie2n0=[0;0;0]; % earth rotation
    case 3
        % test on inertial data captured by the microstrain IMU mounted on a car.
        R_H2_Im = [0, 1, 0; -1, 0, 0; 0, 0, 1];
        p_H2_Im = [0; 0; -(123.6 - 102.8) - 17.8 / 2 - 3] * 1e-2; % +/- 5cm 
        p_H1_Af = [1.15; 0.48; -1.23]; % +/- 2cm
        p_H1_Ar = [-0.75; 0.45; -1.20]; % +/- 2cm
        R_H2_H1 = eye(3);
        p_H2_H1 = [0; -0.21; 0]; % +/- 2cm

        isOutNED=true;
        resdir='/home/jhuai/Desktop/temp/gnssimu/';
        filresfile=[resdir, 'filresult.csv']; % navigation states
        imuresfile=[resdir, 'imuresult.csv']; % imu error terms
        % imu options
        options.startTime=327600;
        options.endTime=options.startTime + 400; % 338555;

        options.imuErrorModel=5; % 4 for random constant acc bias and gyro bias, 5 for random walk acc bias and random constant gyro bias
        options.mechanization=2;

        options.imutype=5; % 5 for microstrain
        options.dt = 0.06;
        options.maxCovStep=1/30; %maximum covariance propagation step, if equal to dt, means single speed mode

        imufile='/media/jhuai/SeagateData/jhuai/data/osu-spin-lab/MultiSensor_NOV_11_2015/IMU_MicroStrain/3DM-GX3 Data Log D.csv';

        imuDataReader = MicrostrainImuDataReader(imufile, options.startTime);

        % position of the rear antenna at startTime
        inixyz_ant=[592575.6758  -4856608.0160   4078416.7701]';
        options.inillh_ant=ecef2geo_v000(inixyz_ant,0);
        options.Vn=[0;0;0]; % velocity in the navigation frame at startTime
        options.Cb2imu=[0 -1 0; 1 0 0; 0 0 1]; % vehicle body frame to Microstrain IMU frame
        options.qb2n= roteu2qr('xyz', [0.008056641; -0.00378418; 8.948364 * pi / 180]); % using H764G-2 INS_ROLL, INS_PITCH, INS_AZ

        options.InvalidateIMUerrors=false;
        options.initAttVar=2*pi/180; % 1 deg std for roll and pitch, 2 times 1 deg for yaw std
        
        % GPS options
        useGPS=true;
        useGPSstd=true; % use the std in the rtklib GPS solutons
        options.imu_p_ant = R_H2_Im' * (p_H2_H1 - p_H2_Im + R_H2_H1 * p_H1_Ar);
        % gps start and end time
        gpsSE=options.startTime+[0, 350];

        gpsfile='/media/jhuai/SeagateData/jhuai/data/osu-spin-lab/MultiSensor_NOV_11_2015/GPS_rover_solution_best/Rear_antenna.pos';
        allGpsData=loadAllGPSData(gpsfile, [options.startTime, options.endTime], 'lla'); 
%         gpsDataReader = GpsDataReader(gpsfile, options.startTime, 'lla');
        isConstantVel=false; % use constant velocity model
        options.velNoiseStd=1; % velocity noise density m/s^2 in a horizontally axis
        
        % ZUPT options
        % ZUPT start and end Time, n x 2 matrix, n is the number of segment        
        zuptSE=[options.startTime, options.startTime+30];
        options.sigmaZUPT = 0.1;% unit m/s
        numImuDataPerZUPT=round(sqrt(1/options.dt));
        % NHC options
        options.sigmaNHC = 0.1; % unit m/s
        numImuDataPerNHC = inf; % Turn off NHC.
        % minimum velocity before applying the NHC, this option decouples ZUPT and NHC
        options.minNHCVel= 2.0;

        gn0= [0; 0; 9.8008]; % gravity in NED frame at start
        wie2n0=[0;0;0]; % earth rotation

        % state initialization options. 
        % Biases computed by averaging leads to slightly worse results
        % thatn zeros, maybe because the IMU measures the car engine
        % vibration at the beginning.
%         staticUntilIndex = find(allImuData(:, 1) < zuptSE(2), 1, 'last' );
%         gyroBias = mean(allImuData(1:staticUntilIndex, 5:7), 1)';  
%         accelBias = mean(allImuData(1:staticUntilIndex, 2:4), 1)' + options.Cb2imu * quatrot_v000(options.qb2n, gn0, 1);  
%         options.imuErrors = [accelBias; gyroBias];

        options.imuErrors=zeros(6,1);
    otherwise
        error(['Unsupported testing case ' num2str(experim) '!']);
end

fprintf('\nEKF of MEMS IMU on data of experiment %d!\n\n', experim);
numUsedGnssData = 0;
lastUsedGnssDataTime = 0;

numUsedNHC = 0;
numUsedZUPT = 0;
filter = EKF_filter_s0frame_bias(options);

imuIndex = 0;
previmudata = imuDataReader.previous();
preimutime = previmudata(1);
imudata = imuDataReader.current();

gpsIndex=1;
imuCountSinceGnss=1;   % to count how many IMU data after the latest GPS observations
% read the GPS data and align the GPS data with the imu data
if useGPS
    gpsdata = allGpsData(gpsIndex, :)';
    assert(gpsdata(1) >= gpsSE(1));
    gnssStartTime = allGpsData(gpsIndex, 1);
else
    gpsdata=inf;
end

% Start the main INS
initime=preimutime;
imuaccum=zeros(6,1);    % record the accumulated imu measurements
curimutime=imudata(1);
covupt_time=preimutime; %the time that we last updated the covariance

fprintf('Estimating trajectory...\n');
ffilres=fopen(filresfile,'Wt'); % navigation states
fimures=fopen(imuresfile,'Wt'); % imu errors

while (curimutime<options.endTime)
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
        formatString = '%.8f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f\n';
        fprintf(fimures, formatString, [curimutime;filter.imuErrors;...
            full(sqrt(diag(filter.p_k_k(filter.imuBiasDriftSIP+(0:5),...
            filter.imuBiasDriftSIP+(0:5)))))]);
    end
    %% ZUPT (zero velocity updates).
    isStatic =~isempty(zuptSE) && ~isempty(find(((zuptSE(:,1)<=curimutime)&(zuptSE(:,2)>=curimutime))==1,1));
    isZUPT = numImuDataPerZUPT < imuCountSinceGnss && mod(imuCountSinceGnss, numImuDataPerZUPT)==1;
    if (isStatic&&isZUPT)
        measure=zeros(3,1);
        predict=filter.rvqs0(4:6);
        H=sparse([zeros(3) eye(3) zeros(3,size(filter.p_k_k,1)-6)]);
        R=eye(3)*options.sigmaZUPT^2;
        filter.correctstates(predict,measure, H, R, curimutime, false);
        numUsedZUPT = numUsedZUPT + 1;
    end

    %% NonHolonomic constraints.
    isNHC = numImuDataPerNHC < imuCountSinceGnss && mod(imuCountSinceGnss, numImuDataPerNHC)==1;
    velnorm=filter.GetVelocityMag();
    if (isNHC&&velnorm>options.minNHCVel)
        % non-holonomic constraints
        Cs02s=quat2dcm_v000(filter.rvqs0(7:10));
        Cs02b=options.Cb2imu'*Cs02s;
        predict=Cs02b(2:3,:)*filter.rvqs0(4:6);
        larry=options.Cb2imu'*skew(quatrot_v000(filter.rvqs0(7:10),filter.rvqs0(4:6),0));
        H=sparse([zeros(2,9) Cs02b(2:3,:) larry(2:3,:) zeros(2,size(filter.p_k_k,1)-15)]);
        
        measure=[0;0];
        R=eye(2)*options.sigmaNHC^2;
        filter.correctstates(predict,measure, H, R, curimutime, false);
        numUsedNHC = numUsedNHC + 1;
    end

    %% GNSS observations.
    if (curimutime>=gpsdata(1))
        if mod(gpsIndex, 20) == 0
            fprintf('Using GNSS data at %.3f.\n', gpsdata(1, 1)); 
        end
        imuCountSinceGnss=0;
        % get antenna position in N frame.
        % 1. matlab mapping toolbox
        [north, east, down] = geodetic2ned(gpsdata(2), gpsdata(3), gpsdata(4), ...
            options.inillh_ant(1), options.inillh_ant(2), options.inillh_ant(3), ...
            wgs84Ellipsoid, 'radians');
        p_N_ant = [north; east; down];
        % 2. our own implementation
%         gpsecef=lla2ecef([gpsdata(2:3) * 180 / pi; gpsdata(4)]')';
%         p_N_ant_goodold = quatrot_v000(filter.rqs02e(4:7), gpsecef - inixyz_ant, 1);
%         assert(max(abs(p_N_ant - p_N_ant_goodold)) < 1e-8);
        % 3. matlab automated driving toolbox.
%         [east2, north2, up2] = latlon2local(gpsdata(2) * 180 / pi, gpsdata(3) * 180 / pi, ...
%             gpsdata(4), [options.inillh_ant(1:2) * 180 / pi; options.inillh_ant(3)]);
%         assert(max(abs(p_N_ant - [north2; east2; -up2])) < 1e-8);

        measure = p_N_ant - quatrot_v000(filter.rvqs0(7:10),options.imu_p_ant,1);

        predict=filter.rvqs0(1:3);

        H=sparse([eye(3), zeros(3), zeros(3,size(filter.p_k_k,1)-6)]);
        if(~useGPSstd)
            R=defaultGpsCovariance(gpsdata(5));
        else
            R=diag((4 * gpsdata(7:9)).^2);
        end
        if gpsdata(1) - lastUsedGnssDataTime > 5 || gpsdata(1) - gnssStartTime < 10
            gatingtest = false;
        end
        applied = filter.correctstates(predict,measure, H, R, gpsdata(1), gatingtest);
        if applied
            lastUsedGnssDataTime = gpsdata(1);
        end

        numUsedGnssData = numUsedGnssData + 1;
        %Read the next gps data that is within the specified sessions
        if gpsIndex < size(allGpsData, 1)
            gpsIndex = gpsIndex + 1;
            gpsdata = allGpsData(gpsIndex, :)';
            if gpsdata(1) >= gpsSE(2) || gpsIndex == size(allGpsData, 1)
                disp(['GNSS outage starts from ' num2str(gpsdata(1)) ...
                    ' GTOW sec which is ', num2str(gpsdata(1) - options.startTime), ...
                    ' since the start!']);
                gpsIndex = size(allGpsData, 1);
                gpsdata(1) = inf;
            end
        end
    else
        imuCountSinceGnss=imuCountSinceGnss+1;
    end
    preimutime=curimutime;
    imuIndex = imuIndex + 1;
    imudata = imuDataReader.next();
    curimutime = imudata(1);
end

fprintf(['EKF filter uses %d IMU data, %d ZUPT constraints, ' ...
    '%d NHC constraints, and %d GNSS observations.\n'], ...
    imuIndex, numUsedZUPT, numUsedNHC, numUsedGnssData);
fclose all;

%% Plot navigation results
kf = readmatrix(filresfile);
posdata=allGpsData;
for i=1:size(posdata,1)
    [north, east, down] = geodetic2ned(posdata(i, 2), posdata(i, 3), posdata(i, 4), ...
            options.inillh_ant(1), options.inillh_ant(2), options.inillh_ant(3), ...
            wgs84Ellipsoid, 'radians');
    posdata(i, 2:4) = [north, east, down];
end

nextFig=1;
f(nextFig) = figure;
plot3(kf(:,1)-kf(1,1),kf(:,3),kf(:,2),'g.');
hold on;
plot3(posdata(:,1)-kf(1,1),posdata(:,3),posdata(:,2),'r+');
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
plot3(posdata(:,2), posdata(:,3), posdata(:, 4), 'r+');
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
plot(posdata(:,1)-kf(1,1),posdata(:,4),'r+');

grid
xlabel('Time [s]')
ylabel('Height/m')
title('Height of antenna by KF(green) and reference(red)');
saveas(f(nextFig),[resdir 'red truth and height'],'fig');

plotkf_v001(kf, resdir);
err = readmatrix(imuresfile);

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

nextFig=nextFig+1;
nextDim=nextDim+3;
f(nextFig) = figure;
loc = nextDim+(0:2);
plot(err(:,1)-err(1,1),err(:,loc),'marker','.')
grid
xlabel('Time [s]')
legend('X','Y','Z')
title('Accelerometer bias drift cov std');

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

