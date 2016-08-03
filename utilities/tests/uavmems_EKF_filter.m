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

function uavmems_EKF_filter()
workspace_path = 'C:\JianzhuHuai\GPS_IMU\programs\matlab_ws';
addpath([workspace_path '\instk']); % imu functions
addpath([workspace_path '\voicebox']); % rotation functions
addpath([workspace_path '\ekfmonocularslamv02']); % filter classes
addpath([workspace_path '\utilities']); % data readers

% than dcm2quat_v000 and dcm2euler_v000
import java.util.LinkedList
clear variables;
clc; close all; format longg;
fprintf('\n EKF of MEMS IMU on UAV!\n\n');
rng('default');

experim=6;
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
        gpspostype=3;           % 1 for calender time 2 for GPSTOW format, both produced by RTKlib
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
        % test on Microstrain 3dm gx3-35 data integration with GPS data
        % collected by octoptor on 2014/10/08 Nikon test       
        isOutNED=true;
        resdir='C:\Users\huai.3\Desktop\huai work\OctoptorINSGPStest\temp\';
        filresfile=[resdir, 'filresult.bin']; % navigation states
        imuresfile=[resdir, 'imuresult.bin']; % imu error terms
        % imu options
        options.startTime= 324310.047; % for Nikon test just after launching which was visually estimated from Nikon images
%          options.startTime= 320335.7513; % for velodyne test just after launching
        options.endTime=323946.4722+262.7+294; % for Nikon test excluding landing
%         320510; % exclude landing for velodyne test
        options.imuErrorModel=5; % 4 for random constant acc bias and gyro bias, 5 for random walk acc bias and random constant gyro bias
        options.mechanization=2; % 1 for wander azimuth
        
        options.imutype=5; % 5 for 3dm gx3 35, 7 for epson m-g362pdc1
        options.dt=1/100;  %sampling interval
        options.maxCovStep=options.dt*2.5; %maximum covariance propagation step, if equal to dt, means single speed mode
        options.imufile='G:\3_Nikon\IMU_MicroStrain-Nikon.csv';
%         'F:\oct082014velodyne\IMU_MicroStrain-Velodyne.csv';
        
        imuFileType=4; % 4 for 3DM GX3 -35 csv, 5 for m-g362pdc1 csv
        inixyz_ant=[579144.429900167         -4851437.31394688          4086499.41041637]'; % for Nikon
%         [ 579144.6076            -4851437.6245           4086499.4477]';% the position of antenna at startTime for velodyne
        options.inillh_ant=ecef2geo_v000(inixyz_ant,0);
        Ce2n0=llh2dcm_v000(options.inillh_ant(1:2),[0;1]);
        
        options.imuErrors=zeros(6,1);
%         options.imuErrors(4:6)=[0.000328736153259098      0.000549543824721639     -0.000343360505567371]'; % for microstrain
%         options.imuErrors(4:6)=[0.00135256381987421       0.00194848735437637      0.000783913812673864]'; % for epson
        options.imuErrors(4:6)=[0.000179418067698059       0.00043635486969162     -0.000164707153981932]'; % for microstrain nikon
%         options.imuErrors(4:6)=[ 6.17841226894643e-05      0.000830677715804075     -0.000116766202634011]'; % for microstrain velodyne
        options.imuErrors(3)=1.5;
%         options.imuErrors(4:6)=[  0.0011371102232815       0.00165615518054893      0.000756598052133374]'; % for epson velodyne
        
        options.Vn=[0;0;0]; % roughly estimated from GPS
        options.qb2n=  roteu2qr('xyz',[0;0; 80]/180*pi); % the rotation from imu to w-frame(rotated n-frame)
%         dcm2quat_v000([0,1,0;1,0,0;0,0,-1]); % from epson frame to microstrain frame
%         dcm2quat_v000(euler2dcm_v000([0.01072662;	-0.00419594; -0.76]));
                
        options.InvalidateIMUerrors=false;
        options.initAttVar=2*pi/180; % 1 deg std for roll and pitch, 2 times 1 deg for yaw std
        
        % GPS options
        useGPS=true; % if to only use height observation, modify if(1) to if(0) in the following code
        useGPSstd=true; % use the std in the rtklib GPS solutons
        % for velodyne microstrain, useGPSstd gives better integration results
        Tant2imu=[8;0;-118]*1e-3; % from Greg's drawing, the position of antenna in microstrain frame
%         Tant2imu=[4;5;99]*1e-3; % from Greg's drawing, the position of antenna in epson frame
        % gps start and end time
        gpsSE=[options.startTime, options.startTime+2000];
        gpspostype=3;           % 1 for calender time 2 for GPSTOW format, both produced by RTKlib
        gpsfile='G:\3_Nikon\GPS_OEM615-Nikon\Processed_respect_Van_base_station\Nikon.pos';
%         'F:\oct082014velodyne\Velodyne.pos';
        isConstantVel=false; % use constant velocity model
        options.velNoiseStd=1; % velocity noise density m/s^2 in a horizontally axis
        
        % ZUPT start and end Time, n x 2 matrix, n is the number of segment
    
%         zuptSE=[options.startTime, 315341.8429+800; 316365.031277514, 316359.031277514+300]; % for oct08 2014 navigation dataset
        zuptSE=[options.startTime, 323946.4722+359];%  320541.8429,  320688]; % for oct08 2014 Nikon dataset
        options.sigmaZUPT = 0.05;% unit m/s
        rateZUPT=round(sqrt(1/options.dt)); % from tests, we see zupt has adverse effect in this case
        % NHC options
        options.sigmaNHC = 0.1;% unit m/s
        rateNHC=inf; % turn off NHC
        % minimum velocity before applying the NHC, this option decouples ZUPT and NHC
        options.minNHCVel=2.0;
        
       % gravity correction related parameters
        rateGravity=3; 
%         gn0=-[ -0.0626187978334586;        -0.095531558468761;         -9.80835081105972];% estimated from static session, may include accelerometer bias for microstrain
%         gn0=-[  -0.00757319038128657        0.0516210641328184          9.80663694349179]'; % for epson
        gn0=-[-0.104112477122421         0.259826191921763          -9.8041828151544]'; % for microstrain nikon test
%         gn0=-[ -0.180653756542627        0.0836986162201811     -9.80614572647246]'; % for microstrain velodyne test
        
        sigmaGravity=2; % unit m/s^2, gravity measurement noise
        gravityMagCutoff=2; 
        wie2n0=[0;0;0]; % no earth rotation
        % vertical velocity constraints
        rateVel=inf; %round(sqrt(1/options.dt)/2);
        sigmaVertVel=6; % unit m/s
        sigmaHorizVel=10; 
               
    case 3
        % test on Epson data integration with GPS data
        % collected by octoptor on 2014/10/08 velodyne and Nikon test       
        isOutNED=true;
        resdir='C:\Users\huai.3\Desktop\huai work\OctoptorINSGPStest\temp\';
        filresfile=[resdir, 'filresult.bin']; % navigation states
        imuresfile=[resdir, 'imuresult.bin']; % imu error terms
        % imu options
        options.startTime=324310.047; % for Nikon test just after launching
        options.endTime=323946.4722+262.7+294; % for Nikon test excluding landing
%         options.startTime=320335.7513; % for velodyne test just after launching
%         options.endTime=320510; % for velodyne test exclude landing
        options.imuErrorModel=5; % 4 for random constant acc bias and gyro bias, 5 for random walk acc bias and random constant gyro bias
        options.mechanization=2; % 1 for wander azimuth
        
        options.imutype=7; % 5 for 3dm gx3 35, 7 for epson m-g362pdc1
        options.dt=1/125;  %sampling interval
        options.maxCovStep=options.dt; %maximum covariance propagation step, if equal to dt, means single speed mode
%         options.imufile='F:\oct082014velodyne\IMU_Epson-tagged-Velodyne-20141008_125016.csv';
        options.imufile='F:\3_Nikon\IMU_Epson-tagged-Nikon-20141008_135928.csv';
        imuFileType=5; % 4 for 3DM GX3 -35 csv, 5 for m-g362pdc1 csv
%         inixyz_ant=[ 579144.6076            -4851437.6245           4086499.4477]';% the position of antenna at startTime for velodyne test
        inixyz_ant=[579144.429900167         -4851437.31394688          4086499.41041637]'; % for Nikon
        options.inillh_ant=ecef2geo_v000(inixyz_ant,0);
        Ce2n0=llh2dcm_v000(options.inillh_ant(1:2),[0;1]);
        
        options.imuErrors=zeros(6,1);
%         options.imuErrors(4:6)=[0.000328736153259098
%         0.000549543824721639     -0.000343360505567371]'; % for
%         microstrain navigation
%         options.imuErrors(4:6)=[0.00135256381987421
%         0.00194848735437637      0.000783913812673864]'; % for epson
%         navigation
%         options.imuErrors(4:6)=[ 6.17841226894643e-05      0.000830677715804075     -0.000116766202634011]'; % for microstrain velodyne
        options.imuErrors(3)=-1.3;
        options.imuErrors(4:6)=[0.00139318997521938       0.00169263883612614      0.000825222230639221]';% for espon nikon test
%         options.imuErrors(4:6)=[  0.0011371102232815       0.00165615518054893      0.000756598052133374]'; % for epson velodyne
        
        options.Vn=[0;0;0]; % roughly estimated from GPS
        options.qb2n= rotro2qr(R3(-80*pi/180)*[0,1,0;1,0,0;0,0,-1]); % the rotation from imu to n-frame
%         dcm2quat_v000([0,1,0;1,0,0;0,0,-1]); % from epson frame to microstrain frame
%         dcm2quat_v000(euler2dcm_v000([0.01072662;	-0.00419594; -0.76]));
        
        options.InvalidateIMUerrors=false;
        options.initAttVar=2*pi/180; % 1 deg std for roll and pitch, 2 times 1 deg for yaw std
        
        % GPS options
        useGPS=false; 
        useGPSstd=true; % use the std in the rtklib GPS solutons
        % for velodyne microstrain, useGPSstd gives better integration results        
        Tant2imu=[4;5;99]*1e-3; % from Greg's drawing, the position of antenna in epson frame
        % gps start and end time
        gpsSE=[options.startTime, options.startTime+2000];
        gpspostype=3;           % 1 for calender time 2 for GPSTOW format, both produced by RTKlib
%         gpsfile='F:\oct082014velodyne\Velodyne.pos';
        gpsfile='F:\3_Nikon\GPS_OEM615-Nikon\Processed_respect_Van_base_station\Nikon.pos';
        
        isConstantVel=false; % use constant velocity model
        options.velNoiseStd=1; % velocity noise density m/s^2 in a horizontally axis
                
        % ZUPT options
        % ZUPT start and end Time, n x 2 matrix, n is the number of segment        
%         zuptSE=[options.startTime, 315341.8429+800; 316365.031277514, 316359.031277514+300]; % for oct08 2014 navigation dataset
%         zuptSE=[options.startTime, 320301.8429];%  320541.8429,  320688]; % for oct08 2014 velodyne dataset
        zuptSE=[options.startTime, 323946.4722+359];%  % for oct08 2014 Nikon dataset
        options.sigmaZUPT = 0.1;% unit m/s
        rateZUPT=round(sqrt(1/options.dt)); % from tests, we see zupt has adverse effect in this case
        % NHC options
        options.sigmaNHC = 0.1;% unit m/s
        rateNHC=inf; % turn off NHC
        % minimum velocity before applying the NHC, this option decouples ZUPT and NHC
        options.minNHCVel=2.0;
        
        % gravity correction related parameters
        rateGravity=1;
%         gn0=-[ -0.0626187978334586;        -0.095531558468761;
%         -9.80835081105972];% estimated from static session, may include
%         accelerometer bias for microstrain navigation test
%         gn0=-[  -0.00757319038128657        0.0516210641328184
%         9.80663694349179]'; % for epson navigation test
%         gn0=quatrot_v000(options.qb2n, -[0.0590959349840854        -0.227331831539279          9.80786160774546]', 0); % for epson velodyne test, gravity in n0 frame
        gn0=quatrot_v000(options.qb2n, -[ 0.237649796494726        -0.155077633109062          9.80415364440748]', 0); % for epson Nikon test, gravity in n0 frame
        sigmaGravity=2; % unit m/s^2, gravity measurement noise
        gravityMagCutoff=2; 
        
        wie2n0=[0;0;0]; % no earth rotation
        
        % vertical velocity constraints
        rateVel=inf; %round(sqrt(1/options.dt)/2);
        sigmaVertVel=6; % unit m/s
        sigmaHorizVel=10; 
        
    case 4
        % test on Epson data integration with GPS data
        % collected by GPS Van on 2015/11/11, session B    
        isOutNED=true;
        resdir='C:\Users\huai.3\Desktop\huai work\OctoptorINSGPStest\temp\';
        filresfile=[resdir, 'filresult.bin']; % navigation states
        imuresfile=[resdir, 'imuresult.bin']; % imu error terms
        % imu options
        options.startTime=321011.2039;
        options.endTime=322027;

        options.imuErrorModel=5; % 4 for random constant acc bias and gyro bias, 5 for random walk acc bias and random constant gyro bias
        options.mechanization=2; % 1 for wander azimuth
        
        options.imutype=7; % 5 for 3dm gx3 35, 7 for epson m-g362pdc1
        options.dt=1/500;  %sampling interval
        options.maxCovStep=5*options.dt; %maximum covariance propagation step, if equal to dt, means single speed mode

        options.imufile='\\File.ceegs.ohio-state.edu\SPIN\MultiSensor_NOV_11_2015\IMU_Epson\20151111_120956_B_timetagged.csv';
        imuFileType=5; % 4 for 3DM GX3 -35 csv, 5 for m-g362pdc1 csv
        % the position of antenna at startTime
        inixyz_ant=[592574.4116  -4856603.9732   4078414.5768]';
        options.inillh_ant=ecef2geo_v000(inixyz_ant,0);
        
        options.imuErrors=zeros(6,1);
        
        options.Vn=[0;0;0]; % velocity in n-frame at start
        options.Cb2imu=[1 0  0;  0 -1 0; 0 0 -1]; % from vehicle body frame FDR to sensor frame
        Cimu2n=[1 , 0, 0; 0, -1,0; 0, 0, -1];
        options.qb2n= dcm2quat_v000(Cimu2n*options.Cb2imu);

        options.InvalidateIMUerrors=false;
        options.initAttVar=2*pi/180; % 1 deg std for roll and pitch, 2 times 1 deg for yaw std
        
        % GPS options
        useGPS=true; 
        useGPSstd=true; % use the std in the rtklib GPS solutons
        Tant2imu=[0.454; -0.746; 1.344];
        % gps start and end time
        gpsSE=options.startTime+[100, 200; 250, 400; 450, 600; 650, 800; 950, 2000];
        gpspostype=4;           % 1 for calender time 2 for GPSTOW format, both produced by RTKlib
        gpsfile='C:\JianzhuHuai\GPS_IMU\programs\matlab_ws\data\Novatel_rover_KIN_20151111.pos';
        
        isConstantVel=false; % use constant velocity model
        options.velNoiseStd=1; % velocity noise density m/s^2 in a horizontally axis
                
        % ZUPT options
        % ZUPT start and end Time, n x 2 matrix, n is the number of segment        
        zuptSE=[options.startTime, options.startTime+100];
        options.sigmaZUPT = 0.1;% unit m/s
        rateZUPT=round(sqrt(1/options.dt));
        % NHC options
        options.sigmaNHC = 0.1;% unit m/s
        rateNHC=inf; % turn off NHC
        % minimum velocity before applying the NHC, this option decouples ZUPT and NHC
        options.minNHCVel=2.0;
        
        % gravity correction related parameters
        rateGravity=10;
        gn0= [ 0; 0; 9.804]; % gravity in the n-frame at start
        sigmaGravity=2; % unit m/s^2, gravity measurement noise
        gravityMagCutoff=2; 
        
        wie2n0=[0;0;0]; % earth rotation in n frame at start        
        % vertical velocity constraints
        rateVel=inf; %round(sqrt(1/options.dt)/2);
        sigmaVertVel=6; % unit m/s
        sigmaHorizVel=10; 
    case 5
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
        gpspostype=4;           % 1 for calender time 2 for GPSTOW format, both produced by RTKlib
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
    case 6
        % test on smartphone inertial data collected by iPhone6S      
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
        options.dt=1/30;  %sampling interval
        options.maxCovStep=options.dt; %maximum covariance propagation step, if equal to dt, means single speed mode

        options.imufile='\\File.ceegs.ohio-state.edu\SPIN\MultiSensor_NOV_11_2015\IMU_Epson\20151111_140059_D_timetagged.csv';
        imuFileType=5; % 4 for 3DM GX3 -35 csv, 5 for m-g362pdc1 csv
        % the position of antenna at startTime
        inixyz_ant=[592574.6611  -4856604.0417   4078414.4645]';
        options.inillh_ant=ecef2geo_v000(inixyz_ant,0);
   
        options.imuErrors=zeros(6,1);        
        options.Vn=[0;0;0]; % roughly estimated from GPS
        options.Cb2imu=[1 0  0;  0 1 0; 0 0 1];% vehicle body frame to imu frame
        accel= [];
        eulers2n=orientbymagnaccel(accel);
        options.qb2n= roteu2qr('xyz', eulers2n); % estimated from accelerometer

        options.InvalidateIMUerrors=false;
        options.initAttVar=2*pi/180; % 1 deg std for roll and pitch, 2 times 1 deg for yaw std
        
        % GPS options
        useGPS=false; 
        useGPSstd=true; % use the std in the rtklib GPS solutons
        Tant2imu=[0.454; -0.746; 1.344];
        % gps start and end time
        gpsSE=options.startTime+[0, 350; 420, 500; 540, 600];
        gpspostype=4;           % 1 for calender time 2 for GPSTOW format, both produced by RTKlib
        gpsfile='C:\JianzhuHuai\GPS_IMU\programs\matlab_ws\data\Novatel_rover_KIN_20151111.pos';
        
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
        error('Unsupported testing case!');
end
Counter.numimurecords=40; % the number of data put in preimudata
% Initialize the model state and covariance of state, process noise and
% measurment noise

filter =EKF_filter_s0frame_bias(options);

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
    imuaccum=filter.ffun_state(imuaccum,[imudata(2:7);preimutime;curimutime; gn0; wie2n0], 0, isConstantVel);
    
    %Update the covariance
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
    %Apply ZUPT. Don't apply zupt for each imu sample. Zupt
    %for each sample must be performed on nominal trajectory, not with
    %the Kalman filter.
    isStatic =~isempty(zuptSE) && ~isempty(find(((zuptSE(:,1)<=curimutime)&(zuptSE(:,2)>=curimutime))==1,1));
    isZUPT =mod(imuctr, rateZUPT)==0;
    if (isStatic&&isZUPT)
        measure=zeros(3,1);
        predict=filter.rvqs0(4:6);
        H=sparse([zeros(3) eye(3) zeros(3,size(filter.p_k_k,1)-6)]);
        R=eye(3)*options.sigmaZUPT^2;
        filter.correctstates(predict,measure, H,R);
    end
    % apply gravity norm measurement
    useGravity =mod(imuctr, rateGravity)==0;
    if (useGravity&& abs(norm(imudata(2:4))-9.81)< gravityMagCutoff)
        predAccel=imudata(2:4)-filter.imuErrors(1:3);
        predict=-quatrot_v000(filter.rvqs0(7:10),predAccel,1);
        H=[zeros(3,6) rotqr2ro(filter.rvqs0(7:10))'*skew(predAccel) zeros(3,filter.covDim-9)];
        measure=gn0;
        R=sigmaGravity^2*eye(3);
        filter.correctstates(predict,measure, H,R);
    end
    isNHC =mod(imuctr, rateNHC)==0;
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
    %%Apply vertical velocity observation
    useVel =mod(imuctr, rateVel)==0;
    if (useVel)        
        predict=filter.rvqs0(4:6);
        H=[eye(3), zeros(3,filter.covDim-3)];
        measure=zeros(3,1);
        R=diag([sigmaHorizVel,sigmaHorizVel,sigmaVertVel].^2);
        filter.correctstates(predict,measure, H,R, 1);
    end
    
    %%Apply the  GPS observations
    
    if (curimutime>=gpsdata(1))
        imuctr=0; % to count how many imu epochs after the recent gps observations
        gpsecef=gpsdata(2:4);
        
        measure= quatrot_v000(filter.rqs02e(4:7), gpsecef - inixyz_ant, 1)- quatrot_v000(filter.rvqs0(7:10),Tant2imu,1);
        predict=filter.rvqs0(1:3);
        %% use three channels
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
        else %% only use height
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
        [fgps, gpsdata]=grabnextgpsdata(fgps, gpspostype);        
        if (gpsdata(1)>gpsSE(gpsSErow,2))
            disp(['GPS outage starts from ' num2str(lastgpstime) 'GTOW sec!']);
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
                    disp(['GPS completely gone from ' num2str(lastgpstime) 'GTOW sec!']);
                    gpsdata(1)=inf;% completely stop using GPS
                else
                    while(gpsdata(1)< gpsSE(nextSErow,1))
                        [fgps, gpsdata]=grabnextgpsdata(fgps, gpspostype);                        
                        while(gpsdata(1) > gpsSE(nextSErow,2) && nextSErow< size(gpsSE,1))
                            nextSErow= nextSErow+1;                            
                        end
                    end
                    disp(['GPS will resume from ' num2str(gpsdata(1)) 'GTOW sec!']);            
                end
            else
                nextSErow = gpsSErow;
                disp(['GPS resumes from ' num2str(gpsdata(1)) 'GTOW sec!']);
            end         
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
kf = readdata(filresfile, 1+18);

% Extracting GPS data
posdata=loadAllGPSData(gpsfile, [options.startTime, options.endTime], gpspostype); % reads specified columns
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

plot3(kf(:,1)-kf(1,1),kf(:,3),kf(:,2),'g.')
hold on 
plot3(posdata(:,1)-kf(1,1),xs0(:,2),xs0(:,1),'r+')
grid
axis equal
    xlabel('Time [s]')
    ylabel('East [m]')
    zlabel('North[m]')

title('green KF trajectory and the red GPS reference for GPS antenna');
% saveas(f(nextFig),[resdir,'red truth and track'],'fig');

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

