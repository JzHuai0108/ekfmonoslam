% General state space model for gps imu integration
% The following state space model is used :
% the complete state is defined as delta(ecef x,y,z), delta(vx, vy, vz),
% delta(roll, pitch, yaw) of the imu, accelerometer and gyro bias drift, acc
% scale and gyro scale factor error
% accelerometer and gyro bias modeled as random walk,
% acc scale and gyro scale factor modeled as random walk

% attitude expressed in quaternion, taking 4 dimensions, but in the P
% matrix, their covariance are in terms of radians^2,
% so attitude takes only 3 dimensions. The P matrix is in full form, not
% sparse.

% The state dynamics are driven by a 18 (accel and gyro white noise,
% bias white noise and scale factor white noise) dimensional white Gaussian
% noise and the observations are corrupted by white Gaussian noise
classdef EKF_filter_eframe < handle
    properties (Hidden)
        type = 'ekf';
        tag  = 'EKF_IMU_GPS_EFRM';  % ID tag       
        covDim=9+12; % the dimension of covariance matrix
    end
    % The following properties can be set only by class methods
    properties (SetAccess = private)
        invalidateIMUerrors; % reset IMU errors when they goes beyond expectation
        imuType; % the type of IMU, e.g., MEMS 3DX GM 3-35, HG1700, etc.   
        imuErrorModel=3;
        % the IMU error model corresponding to random walk bias and
        % random walk scale factor error, 
        % 1 without turn on bias estimates,first order GM bias and random walk scale factor errors
        % 2 with random constant turn on bias estimates, first order GM bias and random walk scale factor errors
        % 3 assumes bias and scale factors are random walks
        % 4 assumes random constant bias and random constant scale factor errors
        
        % the body frame, in our case is assumed to be aligned with the H764G frame,
        % but here we are testing on an arbitray sensor, like a mems imu, so Cb2imu
        % represent the rotations between the body frame and the mems imu frame
        % set it as [] to be determined using IMU coarse alignment
        Cb2imu;
        dt; % sampling interval of IMU, unit sec
        % the translation from antenna to mems frame, i.e., the antenna's
        % position in the mems frame. 
        Tant2imu;        
        rvqs2e; % position of IMU in e frame, v in e, and qs2e
        imuErrors; %ba, bg, sa, sg
        
        imuOrientSIP=7; % imu orientation covariance start index
        imuBiasDriftSIP=10; %imu bias error covariance start index
        imuScaleFactorSIP=16;
        p_k_k; % the covariance of the entire state vector
        camPose; % used by KLT tracker
    end
    methods
        function filter = EKF_filter_eframe(options)
            filter.invalidateIMUerrors=options.InvalidateIMUerrors;
            filter.imuType=options.imutype;
            filter.imuErrorModel=options.imuErrorModel;      
            filter.dt=options.dt;
            filter.Cb2imu=eye(3);       
            filter.Tant2imu=filter.Cb2imu*(options.Tant2body-options.Timu2body);   
            filter.camPose=zeros(7,1);
            filter.camPose(5:7)=options.Cimu2cam*filter.Cb2imu*(options.Timu2body-options.Tcam2body);
            filter.camPose(1:4)=rotro2qr(options.Cimu2cam);  

        % initialize states and covariance for imu gps integration,
        % For imu, rs in e, vs in e, q s2e, ba, bg, sa, sg. (s denotes imu)
            % the position of GPS antenna
            inillh=options.inillh_ant;
            Ve=options.Ve;
            % Cb2n initial value given by external reference
            initqbn=options.qb2n;   % H764G INS data, this b refers to the vehicle body frame
            Cen=llh2dcm_v000(inillh(1:2),[0;1]);
            if(isempty(filter.Cb2imu))
                qs2b=getqimu2body(options.imufile, options.inillh_ant, options.startTime, filter.dt);
                filter.Cb2imu=quat2dcm_v000(qs2b)';
                qs2e=quatmult_v001(rotro2qr(Cen),quatmult_v001(initqbn,qs2b,0),1);
            else
                qs2e=rotro2qr(Cen'*quat2dcm_v000(initqbn)*filter.Cb2imu');
            end
            xyz_imu=ecef2geo_v000(inillh,1)-quatrot_v000(qs2e,filter.Tant2imu,0);
            
            % Initial state vector and covariance matrix
            filter.rvqs2e= [xyz_imu; Ve; qs2e];
            filter.imuErrors=options.imuErrors;
            
            filter.p_k_k=zeros(filter.covDim);
            filter.p_k_k(1:3,1:3)=eye(3)*5^2; %position errors in meter
            filter.p_k_k(4:6,4:6)=eye(3)*0.1^2; %vel error in m/s
            filter.p_k_k(filter.imuOrientSIP+(0:2),filter.imuOrientSIP+(0:2))=...
                diag([options.initAttVar,options.initAttVar,options.initAttVar*2].^2);
            
            IMU_ERRDEF=imu_err_defs_v000(options.imutype);
            filter.p_k_k(filter.imuBiasDriftSIP+(0:2),filter.imuBiasDriftSIP+(0:2))=4*eye(3)*IMU_ERRDEF.acc_bias_var;% enlarge initial std by 2
            filter.p_k_k(filter.imuBiasDriftSIP+(3:5),filter.imuBiasDriftSIP+(3:5))=4*eye(3)*IMU_ERRDEF.gyro_bias_var;
            
            filter.p_k_k(filter.imuScaleFactorSIP+(0:2),filter.imuScaleFactorSIP+(0:2))=4*eye(3)*IMU_ERRDEF.acc_scale_var;
            filter.p_k_k(filter.imuScaleFactorSIP+(3:5),filter.imuScaleFactorSIP+(3:5))=4*eye(3)*IMU_ERRDEF.gyro_scale_var;
            
        end
        %===============================================================================================
        %-- State transition function
        % propogate state with accelerometer and gyro input at time k-1 to state at k,
        % i.e., X(k|k-1), from state at k-1. U1 contains IMU measurement
        % acc, and gyro angular rate
        % and 7th row is the previous epoch(k-1) and 8th row is the current epoch(k)
        % imuaccum records the accumulated delta v and delta angle for two speed
        % covariance update
        function imuaccum= ffun_state(filter, imuaccum, U1)
            if(sum(abs(U1(1:3,1)))<1)
                dt1=U1(8,end)-U1(7,1);
                gyroinc=sum(U1(4:6,:),2);
                accinc=sum(U1(1:3,:),2);
                %reset the imu errors if estimates diverges
                IMU_ERRDEF=imu_err_defs_v000(filter.imuType);
                
                if(filter.invalidateIMUerrors)
                    spuracc=find(abs(filter.imuErrors(1:3))>2*IMU_ERRDEF.initacc_bias_err,1);
                    spurgyro=find(abs(filter.imuErrors(4:6))>2*IMU_ERRDEF.initgyro_bias_err,1);
                    if(~isempty(spuracc)||~isempty(spurgyro))
                        disp(['IMU errors diverge at ' num2str(U1(8,end)) '!']);
                        filter.imuErrors=zeros(12,1);
                    end
                end
                
                angleinc=gyroinc-filter.imuErrors(4:6)*dt1-diag(gyroinc)*filter.imuErrors(6+(4:6))/1000;
                velinc=accinc-filter.imuErrors(1:3)*dt1-diag(accinc)*filter.imuErrors(6+(1:3))/1000;
                %imu data accumulator for the covariance update
                imuaccum=imuaccum+[velinc;angleinc];
                gyroinc=angleinc/dt1;
                accinc=velinc/dt1;
            else
                dt1=U1(8,end)-U1(7,1);
                gyroinc=mean(U1(4:6,:),2);
                accinc=mean(U1(1:3,:),2);            
                %reset the imu errors if estimates diverges
                IMU_ERRDEF=imu_err_defs_v000(filter.imuType);
                
                if(filter.invalidateIMUerrors)
                    spuracc=find(abs(filter.imuErrors(1:3))>2*IMU_ERRDEF.initacc_bias_err,1);
                    spurgyro=find(abs(filter.imuErrors(4:6))>2*IMU_ERRDEF.initgyro_bias_err,1);
                    if(~isempty(spuracc)||~isempty(spurgyro))
                        disp(['IMU errors diverge at ' num2str(U1(8,end)) '!']);
                        filter.imuErrors=zeros(12,1);
                    end
                end
                gyroinc=gyroinc-filter.imuErrors(4:6)-diag(gyroinc)*filter.imuErrors(6+(4:6))/1000;
                accinc=accinc-filter.imuErrors(1:3)-diag(accinc)*filter.imuErrors(6+(1:3))/1000;
                %imu data accumulator for the covariance update
                imuaccum=imuaccum+[accinc;gyroinc]*dt1;                
            end
            
            %run strapdown with angle and velocity increments
            [qs2e_new, Ve_new, ecef_new]=strapdown_ecef_quat_v001(filter.rvqs2e(7:10), filter.rvqs2e(4:6), filter.rvqs2e(1:3), accinc, gyroinc, dt1);
            filter.rvqs2e=[ecef_new;Ve_new;qs2e_new];
        end
        function ffun_covariance(filter, imuaccum, covupt_time, curimutime )
            %propagate the covariance corresponds to states, rs in e, vs in e, q s2e,
            % ba, bg, sa, sg
            covdt=curimutime-covupt_time;
            Cs2e=quat2dcm_v000(filter.rvqs2e(7:10));
            [STM, Qd]=sys_ecef_dcm_v001(filter.rvqs2e(1:3), Cs2e, imuaccum(1:3)/covdt,imuaccum(4:6)/covdt, covdt,filter.imuType, filter.imuErrorModel);
            filter.p_k_k=STM*filter.p_k_k*STM'+Qd;  % the covariance of the navigation states and imu error terms
            
            % Propagate the imu error states, in general make no difference so usually commented
            filter.imuErrors=STM(filter.imuBiasDriftSIP:end,filter.imuBiasDriftSIP:end)*filter.imuErrors;
        end
        %==============================================================================================
        %update both the covariance and state
        function correctstates(filter, predict,measure, H,R, type)           
            p_km1_k=filter.p_k_k;
            inno= predict-measure;
            %Kalman
            K=p_km1_k*H'/(H*p_km1_k*H'+R);
            deltaX=K*inno;
            % update covariance
            filter.p_k_k=(eye(size(p_km1_k,1))-K*H)*p_km1_k*(eye(size(p_km1_k,1))-K*H)'+K*R*K';
            % compute updated state            
            [filter.rvqs2e(7:10), filter.rvqs2e(4:6),filter.rvqs2e(1:3)]=...
                correctnav_eframe_v000(filter.rvqs2e(7:10), filter.rvqs2e(4:6),filter.rvqs2e(1:3), deltaX(1:filter.imuOrientSIP+2));
            filter.imuErrors = filter.imuErrors + deltaX(filter.imuBiasDriftSIP:end);
        end
        function SaveToFile(filter, inillh_ant, preimutime, ffilres)
            Cen=llh2dcm_v000(ecef2geo_v000(filter.rvqs2e(1:3,1),0),[0;1]);
            vr_let=quat2dcm_v000(filter.rvqs2e(7:10));
            vr_c=rotro2eu('xyz',Cen*vr_let)*180/pi; % Cs2n          
            xyz_ant=filter.rvqs2e(1:3)+quatrot_v000(filter.rvqs2e(7:10),filter.Tant2imu,0);
            if(~isempty(inillh_ant))
                vr_a=posdiff_v001(xyz_ant, inillh_ant);
%                 fwrite(ffilres,[preimutime;vr_a;filter.rvqs2e(4:6);vr_c;sqrt(diag(filter.p_k_k(1:9,...
%                     1:9)))],'double');
                fwrite(ffilres,[preimutime;vr_a;filter.rvqs2e(4:6);rotqr2eu('xyz', filter.rvqs2e(7:10));sqrt(diag(filter.p_k_k(1:9,...
                    1:9)))],'double');
            else
                fwrite(ffilres,[preimutime;ecef2geo_v000(xyz_ant,0);filter.rvqs2e(4:6);vr_c;sqrt(diag(filter.p_k_k(1:9,...
                    1:9)))],'double');
            end
        end
        function mag=GetVelocityMag(filter)
            mag=norm(filter.rvqs2e(4:6),2);
        end
   end
end