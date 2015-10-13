
close all
clear;
fprintf('\n cam_ins_gps.m to test filtering techniques for combining GPS, mems IMU and camera \n\n');
%--- General setup
% addpath('C:\JianzhuHuai\GPS_IMU\programs\ekfmonocularslam');
experim=1;
switch experim
    case 0%H764G
        options.isOutNED=true;
        options.useGPSstd=false;% use the std from gps solution
        options.dt=1/256;  %sampling time
        startTime = 500600;% 499600; %500699 is the start time of movement, if we use the attitude of INS, divergenece occurs.
        endTime =500800.0;
        %the level arm offset given by direct meauring
        options.LA          = [    -0.746      0.454     -1.344];
        TIME.gpsNum=6*5;% 5Hz* x s, the interval having GPS
        options.imutype=4;%H764G
        options.zuptSE=[500600, 500689];% zupt start and end time
        resdir='C:\Users\huai.3\Desktop\source_testdata_specs\20130719\';
        imufile=[resdir,'h764g071913.txt'];
        gpspostype=1; % 1 for calender time 2 for GPSTOW format, both produced by RTKlib
        gpsfile=[resdir,'rtkout\oem615_20130719.pos'];% right now we use this format,
        %but later we need to change to some internal gps format
        %Initial PVA
        inillh=[39.958499381*pi/180  -83.055230911*pi/180   205.2000]'; %500600 position from GPS
        Vn=[0;0;0];
        iniRPY = [-0.006225586	-0.002166748	-90.32959*pi/180];% 500600 from h764g
        qbn=att2qua(iniRPY);
        options.useCam=false; 
        options.useGPS=true;
        options.InvalidateIMUerrors=false;
        options.initAttVar=1*pi/180;
        gpsfiletrack=[resdir 'rtkout\071913_chen_Topcon2_1.pos'];
        options.sigmaZUPT =0.1;
        options.sigmaNHC =0.1;
        options.sigmaCAM = 10;
        %NHC and ZUPT options
        rateZUPT=round(sqrt(1/options.dt));
        rateNHC=round(sqrt(1/options.dt));
        options.minNHCVel=0.5;
        options.Cimu2body=1;%att2Cbn([-2.69729070186719;-17.0370732271539;0]*pi/180);
        options.kmlfilename=[resdir, 'googleearth.kml'];
    case 1% MEMS SENSOR Personal navigator, mems data 2013 aug
        %MEMS 3dm gx3-35
        options.isOutNED=true;
        options.useGPSstd=true;% use the std from gps solution
        options.dt=1/100;  %sampling time
        startTime=500466.00;
        endTime=500900.00;
        %the level arm offset given by direct meauring
        options.LA=[    -0.1      0.10     -0.80];
        TIME.gpsNum=(680-466)*5;% 5Hz* x s, the interval having GPS
        options.imutype=5;% MEMS 3DM GX 3-35
        options.zuptSE=[500466.4, 500530.00];% zupt start and end time
        resdir='C:\Users\huai.3\Desktop\source_testdata_specs\PN201308\';
        imufile=[resdir,'microstrain_201308PNimu.txt'];
        gpspostype=2; % 1 for calender time 2 for GPSTOW format, both produced by RTKlib
        gpsfile=[resdir,'OEM_GPS.pos'];% right now we use this format,
        inillh=[39.999758506*pi/180  -83.013445516*pi/180   195.7907]';
        Vn=[0;0;0];
        options.useCam=false; 
         options.useGPS=true; 
        qbn=[];%att2qua([0	-20/180*pi 0]);% 500387.2 H764G INS data
        options.InvalidateIMUerrors=false;
        options.initAttVar=5*pi/180;
        gpsfiletrack=gpsfile;
        options.sigmaZUPT =0.1;
        options.sigmaNHC =0.1;
        options.sigmaCAM = 0.05;
        %NHC and ZUPT options
        rateZUPT=round(sqrt(1/options.dt));
        rateNHC=round(sqrt(1/options.dt));
        options.minNHCVel=0.2;
        options.Cimu2body=1;%att2Cbn([-2.69729070186719;-17.0370732271539;0]*pi/180);
        options.kmlfilename=[resdir, 'googleearth.kml'];
    case 2% aug 08 2013 data MEMS 3DM GX 3 35, using gps solution by rtklib by Zhang Xi
        %MEMS 3dm gx3-35
        options.useCam=false; 
        options.useGPS=true; 
        options.isOutNED=true;
        options.useGPSstd=true;% use the std from gps solution
        options.dt=1/100;  %sampling time
        TIME.maxCovStep=1/200; %maximum covariance propagation step
        options.LA=[    -0.76      0.40     -1.21];  
        options.camLA=[2.129+0.04, 0.987+0.02, -(0.949+.03)];
        %Time configuration
        startTime=415500.00; %static since then
        endTime=415275+800.00;% 415275 start moving,  415527.04s having camera measurements
        TIME.gpsNum=100*5;% 5Hz* x s, the interval having GPS
        options.imutype=5;% MEMS 3DM GX 3-35
        options.camtype=3;% Nikon D800
        options.zuptSE=[414840, 415275-100];% zupt start and end time
        resdir='C:\Users\huai.3\Desktop\source_testdata_specs\20130808\';
        imufile=[resdir 'microstrain_20130808imu.txt'];
        gpspostype=2; % 1 for calender time 2 for GPSTOW format, both produced by RTKlib
        gpsfile=[resdir,'oem615_20130809.pos'];
        ptpairfile=[resdir,'nikon080813ptpairs.txt'];
        %Initial PVA, note this is the GPS antenna position, not IMU
        %position
%         inillh=[40.003224879*pi/180  -83.042990034*pi/180  212.6744]';% 414840
        inillh=[40.003784368*pi/180 -83.042904239*pi/180   212.5730]';%415500
        Vn=[0;0;0];
%         qbnprime=att2qua([0.007965088	-0.003631592 5.39978]/180*pi);% H764G INS data, 414840
        qbnprime=att2qua([0.004333496	0.007080078 -80.36499]/180*pi);% H764G INS data 415500, 
        qbn=[];
        options.InvalidateIMUerrors=false;
        options.initAttVar=5*pi/180;
        gpsfiletrack=gpsfile;
        options.sigmaZUPT = 0.1;
        options.sigmaNHC = 0.1;
        options.sigmaCAM =0.0005;
        options.maxPairs=50;
        options.timeshift=0;% time difference
        %NHC and ZUPT options
        rateZUPT=round(sqrt(1/options.dt));
        rateNHC=round(sqrt(1/options.dt));
        options.minNHCVel=.8;
        options.Cimu2body=1;%att2Cbn([-2.69729070186719;-17.0370732271539;0]*pi/180);
        options.kmlfilename=[resdir, 'googleearth.kml'];
        campairs=CQueue();% record the point pairs
   otherwise return;
end

TIME.alignEpoch=40;
preimudata=CQueue();% record the previous imu data
%open files
fimu=fopen(imufile,'r');
fgetl(fimu);% remove the header
hstream= fgetl(fimu);
mass=textscan(hstream,'%f','delimiter',',');
imudata=mass{1};
imudata(2:4,1)=options.Cimu2body*imudata(2:4,1);
imudata(5:7,1)=options.Cimu2body*imudata(5:7,1);
while(imudata(7,1)==0||imudata(1,1)<startTime)
    preimudata.push(imudata(:,end));%record previous imudata
    if(preimudata.size()>TIME.alignEpoch)
        preimudata.pop();
    end
    hstream= fgetl(fimu);
    mass=textscan(hstream,'%f','delimiter',',');
    imudata=mass{1};
    imudata(2:4,1)=options.Cimu2body*imudata(2:4,1);
    imudata(5:7,1)=options.Cimu2body*imudata(5:7,1);
end

if(isempty(qbn))
    databuf=reshape(cell2mat(preimudata.content()),7,preimudata.size());
    %the method of Yuksel does not work out as well as Yudan's coarse
    %alignment
    %     meandat=mean(databuf(2:end,:),2)/options.dt;
    %     %Compute the alignment between the geodetic and the body frame
    %     [Fc wen_n wie_n g]=geoparam_v001(2, Cen(:,3), height, zeros(3,1));
    %     [Cnb E err]=align_opt_v000([meandat(4:6), meandat(1:3)],[[0;0;-g], wie_n] ,[1 0]); %Note:E defined for errors on imu. Therefore, the order of vectors is important
    %     qbn=cbn2quat(Cnb');
    navdata = coarse_alignment(databuf',inillh, [],0, 3);
    options.Cimu2body=att2Cbn([navdata(8:9),0]);% split the rotation vector
    % calibrate the first imudata that is only rotated by identity matrix
    % and to be used for system propagation
    imudata(2:4,1)=options.Cimu2body*imudata(2:4,1);
    imudata(5:7,1)=options.Cimu2body*imudata(5:7,1);
    if(exist('qbnprime','var'))
        qbn=qbnprime;
    else
        databuf([2:4],:)=options.Cimu2body*databuf([2:4],:);
        databuf([5:7],:)=options.Cimu2body*databuf([5:7],:);
        navdata = coarse_alignment(databuf',inillh, [],0, 3);
        iniatt=navdata(8:10)';
        qbn=att2qua(navdata(8:10));% 500387.2 H764G INS data
    end
end

%--- Initialise GSSM

model = gssm_cam_ins_gps('init',options);

% U2=-1, no observation, 1 GPS, 2 ZUPT, 3 NHC
% ftype = input('Inference algorithm  [ srcdkf / pf / sppf / gspf / gmsppf ] : ','s');  %  set type of inference algorithm (estimator) to use :
ftype='srcdkf';
% Initial state
model.inftype=ftype;
Xh =[blh2xyz(inillh)';zeros(3,1);qbn;zeros(6,1)];% 16 dimension,
imugeo=inillh;
Cen=llh2dcm_v000(imugeo(1:2,1),[0,1]);
xyz_ant=Xh(1:3)+Cen'*quatrot_v000(Xh(7:10),options.LA,0);% note the imu position is initialized with GPS antenna,
% but to be consistent that the ini pos diff is 0, we add the LA
inillh_ant=ecef2geo_v000(xyz_ant,0);% antenna position
%=================================================================================================================
%--- Run estimator on observed data (noisy bearing readings)
filresfile=[resdir, 'filresult.bin'];
fimures=fopen([resdir, 'imuresult.bin'],'w');
lastImu=preimudata.back();
ffilres=fopen(filresfile,'w');
if(options.useGPS)
fgps=fopen(gpsfile,'rt');
if(fgps~=-1)
    %%%%%Discard all observation data before the current time
    %remove the header of gps posdata
    hstream= fgetl(fgps);
    while(true)
        if(isempty(strfind(hstream,'%')))
            break;
        else hstream= fgetl(fgps);
        end
    end
    if(gpspostype==1)
        stoic=sscanf(hstream,'%d/%d/%d%d:%d:%f%f%f%f%d%d%f%f%f');
        [weeksec, weeknum]=Calender2GPSWeek(stoic(1:6));
        mess=stoic;
        gpsdata=zeros(1,9); %lattitude longitude in degree and height in meter, Q and no of satels and sdn sde sdu
        gpsdata(1)=weeksec+(mess(4:6)-stoic(4:6))'*[3600;60;1];
        gpsdata(2:9)=mess(7:14);
        
        while (gpsdata(1)<=lastImu(1,end))
            hstream= fgetl(fgps);
            mess=sscanf(hstream,'%d/%d/%d%d:%d:%f%f%f%f%d%d%f%f%f');
            gpsdata(1)=weeksec+(mess(4:6)-stoic(4:6))'*[3600;60;1];
            gpsdata(2:9)=mess(7:14);
        end
    elseif(gpspostype==2)
        stoic=sscanf(hstream,'%d%f%f%f%f%d%d%f%f%f');
        gpsdata=stoic(2:2+9-1); % GPS TOW lattitude longitude in degree and height in meter, Q and no of satels and one reserve
        
        while (gpsdata(1)<=lastImu(1,end))
            hstream= fgetl(fgps);
            stoic=sscanf(hstream,'%d%f%f%f%f%d%d%f%f%f');
            gpsdata=stoic(2:2+9-1);
        end
    end
end
else gpsdata=inf;
end
%load the camera point pairs data
if(options.useCam)
fpair=fopen(ptpairfile,'rt');
if(fpair~=-1)
    %%%%%Discard all observation data before the current time
    %remove the header of point pair data
    hstream= fgetl(fpair);
    while(true)
        if(isempty(strfind(hstream,'%')))
            break;
        else hstream= fgetl(fpair);
        end
    end
    %GPSTOW frm1,GPSTOW frm2, frm1 id,frm2 id, pts1 id, U,V in frm1, pts2 id, U,V of frm2 in pixel(1 based)
    camel=sscanf(hstream,'%f,%f,%d,%d,%d,%f,%f,%d,%f,%f');
    camel(1:2)=camel(1:2)-options.timeshift;
    while (camel(1)<=lastImu(1,end))% we need to ensure that frame 1 has position and attitude estimate
        frm1time=camel(1);
        while(frm1time==camel(1))
            hstream= fgetl(fpair);
            camel=sscanf(hstream,'%f,%f,%d,%d,%d,%f,%f,%d,%f,%f');
            camel(1:2)=camel(1:2)-options.timeshift;
        end
    end
    lastState=[camel(1);zeros(size(Xh))];% to record the state associated with frame1, the last frame,
end
end

% note it is 1 dimension higher than Xh because of timetag
%%%%%Start the main INS
%prepare for the loop
initime=lastImu(1, end);
gpsctr=0;

%256Hz.Therefore, I am going to read 1 records at a time
curimutime=imudata(1,end);
preimutime=lastImu(1, end);

kkk=1;% to count how many imu epochs after the recent gps observations

fprintf('Estimating trajectory...');

        h=model.h;
        hh = h^2;
        %Xh contains the latest state in ecef xyz, vned, quaternion, acc and gyro bias
        Sx = model.Sx;
        % Get and calculate CDKF scaling parameters and sigma point weights
        Xdim  = model.statedim-1;                                % extract state dimension
        Vdim  = model.Vdim;                                    % extract process noise dimension
        Ndim  = model.Ndim;                                    % extract observation noise dimension
        W1 = [(hh - Xdim - Vdim)/hh   1/(2*hh);                  % sigma-point weights set 1
            1/(2*h)                sqrt(hh-1)/(2*hh)];
        nsp1   = 2*(Xdim+Vdim) + 1;          % number of sigma points (first set)
        U2=-1; %no measurement has been received
        
        while (~feof(fimu)&&curimutime<endTime)
            %Write antenna's position and vned, rpy and accel and gyro bias to the files
            
            imugeo=ecef2geo_v000(Xh(1:3),0);
            Cen=llh2dcm_v000(imugeo(1:2,1),[0,1]);
            xyz_ant=Xh(1:3)+Cen'*quatrot_v000(Xh(7:10),options.LA,0);
            vr_c=quat2dcm_v000(Xh(7:10));
            vr_c=dcm2euler_v000(vr_c)*180/pi;
            if(options.isOutNED)
                vr_a=posdiff_v001(xyz_ant', inillh_ant);
                fwrite(ffilres,[preimutime;vr_a;Xh(4:6);vr_c],'double');
            else
                fwrite(ffilres,[preimutime;ecef2geo_v000(xyz_ant,0);Xh(4:6);vr_c],'double');
            end
            
            fwrite(fimures,[preimutime;Xh(11:end)],'double');
            if (imudata(1)-initime)>60
                disp(['Process Time:' num2str(imudata(1))]);
                initime=imudata(1);
            end
            % regenerate the sigma points after measurement or time
            %             update, it is necessary to regenerate the sigma points after
            %             time update as the noise is to be added
            model.Sx=Sx;
            Z=model.generateparticles(model, Xh, []);
            
            U2=-1;% no observation is available yet
            UU1 = cvecrep([imudata(2:7);preimutime;imudata(1)],nsp1);
            
            % TIME UPDATE
            %-- Calculate predicted state mean
            X_ = model.ffun( model, Z(1:Xdim+1,:), Z(Xdim+2:Xdim+Vdim+1,:), UU1);  % propagate sigma-points through process model
            Xh =model.estimation(model, X_, [W1(1,1);W1(1,2)]);
            if(1)% at least this appraoch works very well for time update, though not so rigorous
                tempdiff=model.diffstate(model, X_(:,2:nsp1), cvecrep(X_(:,1),nsp1-1));
                interdiff=tempdiff(:,1:Xdim+Vdim)-tempdiff(:,Xdim+Vdim+1:nsp1-1);
                interdiff(7:9,:)=subangle(tempdiff(7:9,1:Xdim+Vdim),tempdiff(7:9,Xdim+Vdim+1:nsp1-1));
                interadd=tempdiff(:,1:Xdim+Vdim)+tempdiff(:,Xdim+Vdim+1:nsp1-1);
                interadd(7:9,:)=addangle(tempdiff(7:9,1:Xdim+Vdim),tempdiff(7:9,Xdim+Vdim+1:nsp1-1));
            else% the two methods for computing angle difference has no major difference
                tempdiff=X_(:,2:Xdim+Vdim+1)+X_(:,Xdim+Vdim+2:nsp1)-2*cvecrep(X_(:,1),Xdim+Vdim);
                interdiff=model.diffstate(model, X_(:,2:Xdim+Vdim+1),X_(:,Xdim+Vdim+2:nsp1));
                interadd=zeros(size(interdiff));
                
                for hen=1:(Xdim+Vdim)
                    delta1=quatmult_v001(X_(7:10,1+hen),X_(7:10,1),2);
                    delta2=quatmult_v001(X_(7:10,1+Xdim+Vdim+hen),X_(7:10,1),2);
                    interadd(7:9,hen)=quat2rot_v000(quatmult_v001(delta1,delta2,0));% kind of like subangle(qua2att(X1(7:10,i)),qua2att(X2(7:10,i)))
                end
                interadd(1:6,:)=tempdiff(1:6,:);
                interadd(10:end,:)=tempdiff(11:end,:);
            end
            A = W1(2,1) *interdiff;
            B = W1(2,2) *interadd;
            [temp,Sx_] = qr([A B]',0);
            Sx= Sx_';
            
            %Apply ZUPT. Don't apply zupt for each imu sample. Zupt
            %for each sample must be performed on nominal trajectory, not with
            %the Kalman filter.
            isStatic =~isempty(options.zuptSE) && ~isempty(find(options.zuptSE(:,1)<=curimutime && options.zuptSE(:,2)>=curimutime, 1));
            isZUPT =mod(kkk, rateZUPT)==0;
            if (isStatic&&isZUPT&& options.sigmaZUPT>0)
                U2=2; %zupt
                OBS=0;
                W2      = W1;
                W2(1,1) = (hh - Xdim - Ndim(U2))/hh ;                        % sigma-point weights set 2
                [Xh_new,Sx_new]=model.measurementupdate(model,Xh, Sx,OBS,W2,U2);
                Xh=Xh_new;
                Sx=Sx_new;
            end
            
            isNHC =mod(kkk, rateNHC)==0;
            if (isNHC&&options.sigmaNHC>0.0&&sqrt(Xh(4:6)'*Xh(4:6))>options.minNHCVel)% decouple ZUPT and NHC
                % non-holonomic constraints
                U2=3;
                OBS = zeros(2,1);
                W2      = W1;
                W2(1,1) = (hh - Xdim - Ndim(U2))/hh ;                        % sigma-point weights set 2
                [Xh_new,Sx_new]=model.measurementupdate(model,Xh, Sx,OBS,W2,U2);
                Xh=Xh_new;
                Sx=Sx_new;
            end
            %%gps updates
            if (abs(curimutime-gpsdata(1))<options.dt)
                kkk=0; % to count how many imu epochs after the recent gps observations
                OBS=ecef2geo_v000([gpsdata(2)/180*pi;gpsdata(3)/180*pi;gpsdata(4)],1);
                U2=1;
                if(~options.useGPSstd)
                    if(gpsdata(5)==1)
                        model.oNoise{U2}.cov=diag([0.05,0.05,0.15].^2);
                    elseif(gpsdata(5)==2)
                        model.oNoise{U2}.cov=diag([1.0,1.0,2.0].^2);
                    else
                        model.oNoise{U2}.cov=diag([15,15,15].^2);
                    end
                else
                    model.oNoise{U2}.cov=4*diag(gpsdata(7:9).^2);
                end
                W2      = W1;
                W2(1,1) = (hh - Xdim - Ndim(U2))/hh ;                        % sigma-point weights set 2
                [Xh_new,Sx_new]=model.measurementupdate(model,Xh, Sx,OBS,W2,U2);
                Xh=Xh_new;
                Sx=Sx_new;
                %Read the new gps data
                if(~feof(fgps))
                    hstream= fgetl(fgps);
                    if (~ischar(hstream))
                        gpsdata(1)=inf;
                    else
                        if(gpspostype==1)
                            mess=sscanf(hstream,'%d/%d/%d%d:%d:%f%f%f%f%d%d');
                            gpsdata(1)=weeksec+(mess(4:6)-stoic(4:6))'*[3600;60;1];
                            gpsdata(2:6)=mess(7:11);
                        elseif(gpspostype==2)
                            stoic=sscanf(hstream,'%d%f%f%f%f%d%d%f%f%f');
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
            % use point pairs to update
            if(options.useCam)
            if(lastState(2)==0&&abs(curimutime-lastState(1))<options.dt)% initialize the first frame's state
                lastState(2:end)=Xh;
            end
            % update with the second frame matching points, if lastState
            % has been initialized and the current time is approximately
            % frame2 time, % read in the data          
            campairs.empty();            
            if(lastState(2)~=0&&abs(curimutime-camel(2))<options.dt)
                frm1time=camel(1);
                frm2time=camel(2);
                % read in this section until the last line and camel
                % already read in the first line of the next section
                while(frm1time==camel(1))
                    if(campairs.size()<options.maxPairs)
                        campairs.push(camel);
                    end                
                    hstream= fgetl(fpair);
                    if (~ischar(hstream))
                        camel(1)=inf;
                        break;
                    end
                    camel=sscanf(hstream,'%f,%f,%d,%d,%d,%f,%f,%d,%f,%f');
                    camel(1:2)=camel(1:2)-options.timeshift;
                end
                camdata=reshape(cell2mat(campairs.content()),10,campairs.size());
                OBS=zeros(campairs.size(),1);               
                W2      = W1;
                W2(1,1) = (hh - Xdim)/hh;                        % sigma-point weights set 2
                [Xh_new,Sx_new]=model.MU_AdditiveNoise(model,Xh,Sx,OBS,W2,lastState,camdata);
                Xh=Xh_new;
                Sx=Sx_new;
                if(camel(1)-frm2time==0) % no jump
                    lastState=[camel(1);Xh];
                else% jump occurs require new initialization
                    lastState=[camel(1);zeros(size(Xh))];
                end
            end
            end
            %Read the next imu data
            preimudata.push(imudata(:,end));%record previous imudata
            if(preimudata.size()>TIME.alignEpoch)
                preimudata.pop();
            end
            preimutime=curimutime;
            while(preimutime>=curimutime)
                hstream= fgetl(fimu);
                if (~ischar(hstream))
                    imudata=[];
                    break;
                end
                mass=textscan(hstream,'%f','delimiter',',');
                imudata=mass{1};
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
        %---------------------------------------------------------------------------------------------------------

fprintf(' done.\n\n');
fclose all;

kf = readdata(filresfile, 1+9);
if (~options.isOutNED)
    save_google_kml(kf(:,2:4), options.kmlfilename);
end

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
    for i=1:maxl-minl+1
        rovXYZ = blh2xyz(posdata(i,2:4));
        posdata(i,2:4)=posdiff_v001(rovXYZ,inillh_ant);
    end
end

f(21) = figure;
if (options.isOutNED)
%         plot(kf(:,3),kf(:,2),'g.')
%         hold on
%         plot(posdata(:,3),posdata(:,2),'r+')
%         grid
%         xlabel('East [m]')
%         ylabel('North[m]')
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
title('green KF trajectory and the red GPS reference for GPS antenna');
saveas(f(21),[resdir,'red truth and track'],'fig');

f(22) = figure;

plot(kf(:,1)-kf(1,1),kf(:,4),'g.')
hold on
plot(posdata(:,1)-kf(1,1),posdata(:,4),'r+');
grid
xlabel('Time [s]')
ylabel('Height/m')
title('Height of antenna by KF(green) and reference(red)');
saveas(f(22),[resdir 'red truth and height'],'fig')
filresfile=[resdir, 'filresult.bin'];
kf = readdata(filresfile, 1+9);
plotkf(kf, resdir);
err = readdata([resdir, 'imuresult.bin'], 1+6);
ploterr(err,resdir);

