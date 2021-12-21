function forwardfilterMEMS()
%forward filter with internal GPS and IMU data format, primarily designed
%for MEMS sensors

%Each IMU data is arranged as GPS TOW, AX, AY, AZ, WX, WY, WZ
%units sec,m/s^2, m/s^2, m/s^2, rad/s, rad/s,rad/s
% in this program, the imudata is organized as a 7xn matrix, n is the
% number of observations, 7 corresponds to gps time, xyz delta v, xyz delta theta

% this implementation assumes: bias and scale factor residuals are random
% constants. This asuumption makes the model not very suitable for MEMS sensors,
% since they have constant bias, bias drift or instability, scale factor
% instablity, the latter two are better approximated by Gauss Markov first
% order model. This is also indicated by the poor results of this algorithm
% in some cases. For further performance improvement, use better noise
% models.
% The level arm offset is not estimated.

close all;
clear;
caseID=4;
switch caseID          
%     case 1
%         % test on H764G data, to be implemented
%         TIME.dt=1/256;  %sampling time
%         TIME.maxCovStep=1/512; %maximum covariance propagation step
%         options.LA=[ -0.746; 0.454; -1.344];      
%         resdir='C:\Users\huai.3\Desktop\source_testdata_specs\20130719\';
%         imufile=[resdir,'H764GM0719_1.csv'];
%         gpsfile=[resdir,'rtkout\oem615_20130719.pos'];
    case 2    
        resdir='C:\Users\huai.3\Desktop\source_testdata_specs\20130719\';
        options.isOutNED=true;
        options.useGPSstd=false;% use the std from gps solution, otherwise use GPS solution quality based empirical noise covariance,
        options.kmlfilename=[resdir,'0719mems.kml']; 
        TIME.dt=1/100;  %sampling time
        TIME.maxCovStep=1/200; %maximum covariance propagation step
        options.LA=[    -0.76      0.40     -1.21];  
        %Time configuration
        startTime=500387.2;
        endTime=501500;
        imutype=5;% MEMS 3DM GX 3-35
        %--------------------------------------------------------------------------
        % ZUPT start and end Time, determine before the filter
        % format n x 2, n is the number of segment, first  column is the start time
        % second column is the end time, for example [575908      576243; 576502      576602]
        options.zuptSE=[500388, 500600];
        TIME.gpsNum=200*5;% 5Hz* 60s, % 5Hz* x s, the interval having GPS, note 500699 starts movement
        % it seems longer GPS coverage until the car start to move will provide better solution.
        imufile=[resdir,'microstrain_20130719imu.txt'];
        gpspostype=1; % 1 for calender time 2 for GPSTOW format, both produced by RTKlib, but later 2 would be the preferable format
        gpsfile=[resdir,'rtkout\oem615_20130719.pos'];
        options.InvalidateIMUerrors=false;% whether to use UMU correction for MEMS
        options.Cimu2body=1;       
        options.sigmaZUPT = 0.1;    
        % Non-Holonomic Constraint (NHC) sigma (<=0.0 => turn off NHC)
        options.sigmaNHC =0.1; 
        options.minNHCVel=1;% the vehilce moves much faster than human
        %NHC and ZUPT options
rateZUPT=round(sqrt(1/TIME.dt))/2;
rateNHC=round(sqrt(1/TIME.dt))/2;
        inillh=[39.958500453*pi/180;  -83.055237599*pi/180;   205.5195];
        Vn=[0;0;0];
        qbn=att2qua([-0.004486084	-0.002502441 -90.29663*pi/180]);% 500387.2 H764G INS data
        options.initAttVar=2*pi/180;% enlarge this can cause larger drift, bizarre
        gpsfiletrack=[resdir,'rtkout\071913_chen_Topcon2_1.pos'];
    case 3% Personal navigator, mems data 2013 aug
        %MEMS 3dm gx3-35, if we set the bias to 0 after gps outage, the
        %result turns out to be better.
        options.isOutNED=true;
        options.useGPSstd=true;% use the std from gps rtklib output
        TIME.dt=1/100;  %sampling interval
        TIME.maxCovStep=1/200; %maximum covariance propagation step
        options.LA=[    -0.76      0.40     -1.21];  
        %Time configuration
        startTime=500466.00;
        endTime=500900.00;
        TIME.gpsNum=(680-466)*5;% 5Hz* x s, the interval having GPS
        imutype=5;% MEMS 3DM GX 3-35   
        options.zuptSE=[500466.4, 500530.00];% zupt start and end time
        resdir='C:\Users\huai.3\Desktop\source_testdata_specs\PN201308\';
        imufile=[resdir,'microstrain_201308PNimu.txt'];
        gpspostype=2; % 1 for calender time 2 for GPSTOW format, both produced by RTKlib
        gpsfile=[resdir,'OEM_GPS.pos'];
        %Initial PVA
        inillh=[39.999758506*pi/180  -83.013445516*pi/180   195.7907]';
        Vn=[0;0;0];  
        qbn=[];%att2qua([0	-20/180*pi 0]);% 500387.2 H764G INS data
        options.InvalidateIMUerrors=false;
        options.initAttVar=5*pi/180;
        gpsfiletrack=gpsfile;
        options.sigmaZUPT = 0.1;         
        options.sigmaNHC = 0.1;
        %NHC and ZUPT options
rateZUPT=round(sqrt(1/TIME.dt));
rateNHC=round(sqrt(1/TIME.dt));
        options.minNHCVel=0.2;
        options.Cimu2body=1;%att2Cbn([-2.69729070186719;-17.0370732271539;0]*pi/180);
        options.kmlfilename=[resdir, 'googleearth.kml'];
  case 4% aug 08 2013 data MEMS 3DM GX 3 35, using gps solution by rtklib by Zhang Xi
        %MEMS 3dm gx3-35
        options.isOutNED=true;
        options.useGPSstd=true;% use the std from gps solution
        TIME.dt=1/100;  %sampling time
        TIME.maxCovStep=1/200; %maximum covariance propagation step
        options.LA=[    -0.1      0.10     -0.80];  
        %Time configuration
        startTime=415500.00; %static since then
        endTime=415275+800.00;% 415275 start moving,  415527.04s having camera measurements
        TIME.gpsNum=100*5;% 5Hz* x s, the interval having GPS
        imutype=5;% MEMS 3DM GX 3-35   
        options.zuptSE=[414840, 415275-100];% zupt start and end time
        resdir='C:\Users\huai.3\Desktop\source_testdata_specs\20130808\';
        imufile=[resdir, 'microstrain_20130808imu.txt'];
        gpspostype=2; % 1 for calender time 2 for GPSTOW format, both produced by RTKlib
        gpsfile=[resdir,'oem615_20130809.pos'];
        %Initial PVA, note this should be the IMU position but here we use
        % the GPS antenna position to initialize the IMU position, not big
        % deal for kalman filter
%         inillh=[40.003224879*pi/180  -83.042990034*pi/180  212.6744]';
        inillh=[40.003784368*pi/180 -83.042904239*pi/180   212.5730]';%415500
        Vn=[0;0;0];
%         qbnprime=att2qua([0 0 5.39978/180*pi]);% 500387.2 H764G INS data
        qbnprime=att2qua([0.004333496	0.007080078 -80.36499]/180*pi);% H764G INS data 415500, 
        qbn=[];
        options.InvalidateIMUerrors=false;
        options.initAttVar=5*pi/180;
        gpsfiletrack=gpsfile;
        options.sigmaZUPT = 0.1;         
        options.sigmaNHC = 0.1;
        %NHC and ZUPT options
rateZUPT=round(sqrt(1/TIME.dt));
rateNHC=round(sqrt(1/TIME.dt));
        options.minNHCVel=.8;
        options.Cimu2body=1;%att2Cbn([-2.69729070186719;-17.0370732271539;0]*pi/180);
        options.kmlfilename=[resdir, 'googleearth.kml'];
    otherwise
end



TIME.alignEpoch=40;
preimudata=CQueue();% record the previous imu data

%open files
fimu=fopen(imufile,'r');
h= fgetl(fimu);% remove the header
h= fgetl(fimu);
mass=textscan(h,'%f','delimiter',' ');
imudata=mass{1};
imudata(2:7,1)=imudata(2:7,1)*TIME.dt;
imudata(2:4,1)=options.Cimu2body*imudata(2:4,1);
imudata(5:7,1)=options.Cimu2body*imudata(5:7,1);
while(imudata(7,1)==0||imudata(1,1)<startTime) 
    preimudata.push(imudata(:,end));%record previous imudata 
    if(preimudata.size()>TIME.alignEpoch)
        preimudata.pop();
    end
    h= fgetl(fimu);
    mass=textscan(h,'%f','delimiter',' ');
    imudata=mass{1};
imudata(2:7,1)=imudata(2:7,1)*TIME.dt;
imudata(2:4,1)=options.Cimu2body*imudata(2:4,1);
imudata(5:7,1)=options.Cimu2body*imudata(5:7,1);
end
fgps=fopen(gpsfile,'rt');

filresfile=[resdir, 'filresult.bin'];
ffilres=fopen(filresfile,'w');
fimures=fopen([resdir, 'imuresult.bin'],'w');

Cen=llh2dcm_v000(inillh(1:2),[0,1]);
height=inillh(3); %elipsoid height
%Load IMU error definitions
IMU_ERRDEF=imu_err_defs_v000(imutype);
if(isempty(qbn)) 
    if(exist('qbnprime','var'))
        qbn=qbnprime;
    else
        assert(false);
    end
end

%Initial Covariance matrix
nst=21; % Rx RY RZ in earth centered n-frame, vn, RPY, acc bias, gyro bias, acc scale and gyro scale factor residuals
P=zeros(nst);
P(1:3,1:3)=eye(3)*5^2; %position errors in meter
P(4:6,4:6)=eye(3)*0.1^2; %vel error in m/s
P(7:9,7:9)=diag([options.initAttVar,options.initAttVar,options.initAttVar*2].^2); 
P(10:12,10:12)=4*eye(3)*IMU_ERRDEF.acc_bias_var;% enlarge initial std by 2
P(13:15,13:15)=4*eye(3)*IMU_ERRDEF.gyro_bias_var;
P(16:18,16:18)=4*eye(3)*IMU_ERRDEF.acc_scale_var;
P(19:21,19:21)=4*eye(3)*IMU_ERRDEF.gyro_scale_var;

%Discard all observation data before the current time
%remove the header of gps posdata
h= fgetl(fgps); 
while(true)
    if(isempty(strfind(h,'%')))
        break;
    else h= fgetl(fgps);
    end
end
if(gpspostype==1)
    stoic=sscanf(h,'%d/%d/%d%d:%d:%f%f%f%f%d%d%f%f%f');
    [weeksec, weeknum]=Calender2GPSWeek(stoic(1:6));
    mess=stoic;
    gpsdata=zeros(1,9); %lattitude longitude in degree and height in meter, Q and no of satels and sdn sde sdu
    gpsdata(1)=weeksec+(mess(4:6)-stoic(4:6))'*[3600;60;1];
    gpsdata(2:9)=mess(7:14);
    lastImu=preimudata.back();
    while (gpsdata(1)<=lastImu(1,end))
        h= fgetl(fgps);
        mess=sscanf(h,'%d/%d/%d%d:%d:%f%f%f%f%d%d%f%f%f');
        gpsdata(1)=weeksec+(mess(4:6)-stoic(4:6))'*[3600;60;1];
        gpsdata(2:9)=mess(7:14);
    end
elseif(gpspostype==2)
    stoic=sscanf(h,'%d%f%f%f%f%d%d%f%f%f');    
    gpsdata=stoic(2:2+9-1); % GPS TOW lattitude longitude in degree and height in meter, Q and no of satels and one reserve
    lastImu=preimudata.back();
    while (gpsdata(1)<=lastImu(1,end))
        h= fgetl(fgps);
        stoic=sscanf(h,'%d%f%f%f%f%d%d%f%f%f');  
        gpsdata=stoic(2:2+9-1);
    end
end
%%%%%Start the main INS
%prepare for the loop
initime=lastImu(1, end);
covupt_time=lastImu(1, end);  %the time that we last updated the covariance
imuaccum=zeros(6,1);
imuerrors=zeros(12,1);
gpsctr=0;

%read 1 records at a time
curimutime=imudata(1,end);
preimutime=lastImu(1, end);

%Record covariance and navigation solution for the smoother
outimuerrors=[imuerrors(1:3);imuerrors(4:6);imuerrors(7:9)/1000;imuerrors(10:12)/1000];
fwrite(fimures,[curimutime;outimuerrors;sqrt(diag(P(10:21,10:21)))],'double');
kkk=1;% to count how many imu epochs after the recent gps observations
while (~feof(fimu)&&curimutime<endTime)
    %Write result to the files     
    %output GPS antenna position
    vr_c=quat2dcm_v000(qbn); 
    if(options.isOutNED)
        vr_a=posdiff_v000(Cen(3,:), height, inillh)'+quatrot_v000(qbn,options.LA,0);
        fwrite(ffilres,[preimutime;vr_a;Vn;dcm2euler_v000(vr_c)*180/pi;sqrt(diag(P(1:9,1:9)))],'double');
    else 
        xyz_ant=ecef2geo_v000([dcm2llh_v000(Cen);height],1)+Cen'*quatrot_v000(qbn,options.LA,0);  
        fwrite(ffilres,[preimutime;ecef2geo_v000(xyz_ant,0);Vn;dcm2euler_v000(vr_c)*180/pi;sqrt(diag(P(1:9,1:9)))],'double');
    end
    
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
    if(options.InvalidateIMUerrors)
        spuracc=find(abs(imuerrors(1:3))>IMU_ERRDEF.initacc_bias_err,1);      
        spurgyro=find(abs(imuerrors(4:6))>IMU_ERRDEF.initgyro_bias_err,1);    
        if(~isempty(spuracc)||~isempty(spurgyro))
            disp(['IMU errors diverge at ' num2str(curimutime) '!']);
            imuerrors(1:end)=0;
        end   
    end
    angleinc=angleinc-(imuerrors(4:6)*dt+diag(angleinc)*imuerrors(10:12)/1000);
    velinc=velinc-(imuerrors(1:3)*dt+diag(velinc)*imuerrors(7:9)/1000);
    %imu data accumulator for the covariance update
    imuaccum=imuaccum+[angleinc;velinc];
    
    %run strapdown with angle and velocity increments
    [qbn, Vn, Cen, height]=strapdown_Cen_quat_v000(qbn, Vn, Cen, height, velinc, angleinc, dt);
    
    %Update the covariance
    if ((curimutime-covupt_time)>=TIME.maxCovStep) || (curimutime>gpsdata(1))
        %propagate the covariance
        covdt=curimutime-covupt_time;
        % the following two model gives similar results,
        % sys_metric_phipsi_v000 and sys_metric_phipsi_v002
        [STM, Qd]=sys_metric_phipsi_v002(Cen, height, Vn, qbn, imuaccum(4:6)/covdt, imuaccum(1:3)/covdt, covdt,2,imutype); %use psi implementation
      
%        [STM, Qd]=sys_metric_phipsi_v000(Cen, height, Vn, qbn, imuaccum(4:6)/covdt, imuaccum(1:3)/covdt, covdt,2,4); %use psi implementation
%         disp(max(max(abs(STM1-STM))));
%         disp(max(max(abs(Qd1-Qd))));
        P=STM*P*STM'+Qd;
        
        %Propagate the imu error states, make no difference so commented
%         imuerrors=STM(10:end,10:end)*imuerrors; 
        covupt_time=curimutime;
        imuaccum=zeros(6,1);
        
        %Record covariance and navigation solution for the smoother
        %Note that I am recording the predicted solutions. That is why the
        %backward part must use filtered solution.
        outimuerrors=[imuerrors(1:3);imuerrors(4:6);imuerrors(7:9)/1000;imuerrors(10:12)/1000];
        fwrite(fimures,[curimutime;outimuerrors;sqrt(diag(P(10:21,10:21)))],'double');
    end
    %Apply ZUPT. Don't apply zupt for each imu sample. Zupt
    %for each sample must be performed on nominal trajectory, not with
    %the Kalman filter.   
    isStatic =~isempty(options.zuptSE) && ~isempty(find(options.zuptSE(:,1)<=curimutime && options.zuptSE(:,2)>=curimutime, 1));
    isZUPT =mod(kkk, rateZUPT)==0;
    if (isStatic&&isZUPT&& options.sigmaZUPT>0)
        inno=Vn;       
        H=[zeros(3) eye(3) zeros(3,15)];
        R=eye(3)*options.sigmaZUPT^2;        
        %Kalman
        K=P*H'/(H*P*H'+R);
        dx=K*inno;
        P=(eye(21)-K*H)*P*(eye(21)-K*H)'+K*R*K';        
        %Correct states
        [qbn, Vn, Cen, height]=correctnav_Cen_v000(qbn, Vn, Cen, height, dx(1:9), 2, 2);
        imuerrors=imuerrors+dx(10:end);
    end
    
    isNHC =mod(kkk, rateNHC)==0;
    if (isNHC&&options.sigmaNHC>0.0&&sqrt(Vn'*Vn)>options.minNHCVel)% decouple ZUPT and NHC
       % non-holonomic constraints   
        Cnb=quat2dcm_v000(qbn)';
        inno=Cnb(2:3,:)*Vn;      
        H=[zeros(2,3) Cnb(2:3,:) -Cnb(2:3,:)*skew(Vn) zeros(2,12)];
        R=eye(2)*options.sigmaNHC^2;        
        %Kalman
        K=P*H'/(H*P*H'+R);
        dx=K*inno;
        P=(eye(21)-K*H)*P*(eye(21)-K*H)'+K*R*K';
        
        %Correct states
        [qbn, Vn, Cen, height]=correctnav_Cen_v000(qbn, Vn, Cen, height, dx(1:9), 2, 2);
        imuerrors=imuerrors+dx(10:end);    
    end
    
    %Apply the  GPS observations
    if (abs(curimutime-gpsdata(1))<TIME.dt)
        %it is OK to assume spherical earth to form observations (see
        %Savage:15.1.2.2.32)
        %However, If you want to be a little bit more precise use the following
        kkk=0; % to count how many imu epochs after the recent gps observations
        gpsecef=ecef2geo_v000([gpsdata(2)/180*pi;gpsdata(3)/180*pi;gpsdata(4)],1);
        Re=6378137/(sqrt(1.0-0.00669437999014*Cen(3,3)^2));
        posecef=-[(Re+height)*Cen(3,1);(Re+height)*Cen(3,2);(Re*(1-0.00669437999014)+height)*Cen(3,3)];
        inno=Cen*(posecef-gpsecef)+quatrot_v000(qbn,options.LA,0);
        
        H=[eye(3) zeros(3,18)]; %ignore the effect of Cbn errors on level-arm compensation. If you want, use H(:,7:9) = skew(options.LA);
        % the following setting of noise variances is suitable for RTKlib output
        if(~options.useGPSstd)
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
        
        %Kalman
        K=P*H'/(H*P*H'+R);
        dx=K*inno;
        P=(eye(21)-K*H)*P*(eye(21)-K*H)'+K*R*K';

        %Correct states
        [qbn, Vn, Cen, height]=correctnav_Cen_v000(qbn, Vn, Cen, height, dx(1:9), 2, 2);
        imuerrors=imuerrors+dx(10:end);       
        %Read the next gps data 
        if(~feof(fgps))            
            h= fgetl(fgps);
            if (~ischar(h))
                gpsdata(1)=inf;
            else
                if(gpspostype==1)
                    mess=sscanf(h,'%d/%d/%d%d:%d:%f%f%f%f%d%d');
                    gpsdata(1)=weeksec+(mess(4:6)-stoic(4:6))'*[3600;60;1];
                    gpsdata(2:6)=mess(7:11);
                elseif(gpspostype==2)                    
                    stoic=sscanf(h,'%d%f%f%f%f%d%d%f%f%f');
                    gpsdata=stoic(2:2+9-1);
                end
            end
        else
            gpsdata(1)=inf;
        end
        gpsctr=gpsctr+1;
        if (gpsctr>=TIME.gpsNum)
            gpsdata(1)=inf;% stop using GPS, dead reckoning       
        end
    end 
    
    %Read the next imu data
    preimudata.push(imudata(:,end));%record previous imudata 
    if(preimudata.size()>TIME.alignEpoch)
        preimudata.pop();
    end    
    preimutime=curimutime;   
    while(preimutime>=curimutime)
        h= fgetl(fimu);
        if (~ischar(h))
            imudata=[];
            break;
        end
        mass=textscan(h,'%f','delimiter',' ');
        imudata=mass{1};
        imudata(2:7,1)=imudata(2:7,1)*TIME.dt;
        imudata(2:4,1)=options.Cimu2body*imudata(2:4,1);
        imudata(5:7,1)=options.Cimu2body*imudata(5:7,1);
        curimutime=imudata(1,end);
    end
    if (isempty(imudata))
        break;
    else
        kkk=kkk+1;
    end
end

% show the results
fclose all;
kf = readdata(filresfile, 1+9+9);
if (~options.isOutNED)
    save_google_kml(kf(:,2:4), options.kmlfilename);
end
% this section requires that the gps data are of type 2 with GPSTOW
posdata=load(gpsfiletrack);
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
if(options.isOutNED) 
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
end
    
f(21) = figure;
if (options.isOutNED)
    plot(kf(:,3),kf(:,2),'g.')
    hold on
    plot(posdata(:,3),posdata(:,2),'r+')
    grid   
    xlabel('East [m]')
    ylabel('North[m]')
    figure(31)
    plot3(kf(:,1)-kf(1,1),kf(:,3),kf(:,2),'g.')
    hold on
    plot3(posdata(:,1)-kf(1,1),posdata(:,3),posdata(:,2),'r+')
    grid
    axis equal
    xlabel('Time [s]')
    ylabel('East [m]')
    zlabel('North[m]')
else
    plot(kf(:,3),kf(:,2),'g.')
    hold on
    plot(posdata(:,3),posdata(:,2),'r+')
    grid
    axis equal
    xlabel('Lontitude [radian]')
    ylabel('Latitude [radian]')
end
title('Ground Track and red GPS');
saveas(f(21),[resdir,'\red truth and track'],'fig');

f(22) = figure;

plot(kf(:,1)-kf(1,1),kf(:,4),'g.')
hold on
plot(posdata(:,1)-kf(1,1),posdata(:,4),'r+');
grid
xlabel('Time [s]')
ylabel('Height/m')
title('Height of antenna by IMU and red GPS');
saveas(f(22),[resdir 'red truth and height'],'fig')
filresfile=[resdir, 'filresult.bin'];
kf = readdata(filresfile, 1+9+9);
plotkf(kf, resdir);
err = readdata([resdir, 'imuresult.bin'], 1+12+12);
ploterr(err,resdir);
