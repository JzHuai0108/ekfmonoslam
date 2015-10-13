% General state space model for gps imu integration in n-frame
% The following state space model is used :
% the complete state is defined as delta(r x,y,z in earth centered
% "abused" n frame), delta(vx, vy, vz in n frame),
% delta(roll, pitch, yaw) of the imu, accelerometer and gyro bias drift, acc
% scale and gyro scale factor error
% accelerometer and gyro bias are usually modeled as random walk, it can
% also be random constant, This does not make much difference
% acc scale and gyro scale factor are usually modeled as random walk, it can
% also be random constant, This does not make much difference

% attitude expressed in quaternion, taking 4 dimensions, but in the P
% matrix, their covariance are in terms of radians^2,
% so attitude takes only 3 dimensions. The P matrix is in full form, not
% sparse.

% The state dynamics are driven by a 18 (accel and gyro white noise,
% bias white noise and scale factor white noise) or 6 (only accel and
% gyro white noise) dimensional white Gaussian
% noise and the observations are corrupted by white Gaussian noise
% test cases: (1) random walk bias and scale facotr error, random constant
% bias and scale factor error, 
% (2) covariance update step > sample interval
% (3) different imy types, 3DM GX3-35 AND h764g
% (4) psi and phi model
% There is a very subtle point in implementing a class in matlab. Always
% call a class the same name, do not call it a filter in a function and a
% model at another. It causes strange behaviors.
classdef EKF_filter_nframe < handle
    properties (Hidden)
        type = 'ekf';
        tag  = 'EKF_IMU_GPS_NFRM';  % ID tag
        covDim=9+12; % the dimension of covariance matrix, Rx RY RZ in earth centered
        % n-frame, vn, RPY, acc bias, gyro bias, acc scale and gyro scale 
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
        Cen; % rotation from e frame to n frame
        height; % h in geodetic coordinates
        qs2n; % q s(sensor) 2 navigation frame
        Vn; % velocity of the IMU sensor in n frame
        imuErrors; %ba, bg, sa, sg
        
        imuOrientSIP=7; % imu orientation covariance start index
        imuBiasDriftSIP=10; %imu bias error covariance start index
        imuScaleFactorSIP=16;
        p_k_k; % the covariance of the entire state vector
        mode=2; % 1 for phi model, 2 for psi model
        mechanization=2; % 1 for wander azimuth, 2 for local geodetic
        rvqs2e; % added only for compatibility
        camPose; % added only for compatibility
    end
    methods
        function filter = EKF_filter_nframe(options)
            filter.invalidateIMUerrors=options.InvalidateIMUerrors;
            filter.imuType=options.imutype;
            filter.imuErrorModel=options.imuErrorModel;         
            filter.dt=options.dt;
            filter.Cb2imu=options.Cb2imu;
            filter.Tant2imu=filter.Cb2imu*(options.Tant2body-options.Timu2body);
            filter.Vn=options.Vn;         
            % initialize states and covariance for imu gps integration,
            % For imu, rs in n, vs in n, q s2n, ba, bg, sa, sg. (s denotes imu)
            % the position of GPS antenna
            inillh=options.inillh_ant;
            % Cb2n initial value given by external reference
            initqbn=options.qb2n;   % H764G INS data, this b refers to the vehicle body frame
            if(isempty(filter.Cb2imu))
                qs2b=getqimu2body(options.imufile, options.inillh_ant, options.startTime, filter.dt);
                filter.Cb2imu=quat2dcm_v000(qs2b)';
                filter.qs2n=quatmult_v001(initqbn,qs2b,0);
            else
                filter.qs2n=quatmult_v001(initqbn,rotro2qr(filter.Cb2imu),2);
            end
            Ce2n0=llh2dcm_v000(inillh(1:2),[0;1]);
            qs2e=quatmult_v001(rotro2qr(Ce2n0), filter.qs2n,1);
            xyz_imu=ecef2geo_v000(inillh,1)-quatrot_v000(qs2e,filter.Tant2imu,0);
            inillh_imu=ecef2geo_v000(xyz_imu, 0);
            filter.height=inillh_imu(3); %elipsoid height
            filter.Cen=llh2dcm_v000(inillh_imu(1:2),[0;1]);
            filter.imuErrors=options.imuErrors;
            filter.p_k_k=zeros(filter.covDim);
            filter.p_k_k(1:3,1:3)=eye(3)*1^2; %position errors in meter
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
        % delta v, and delta angle of gyro or acc and angular rate
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
                angleinc=gyroinc*dt1;
                velinc=accinc*dt1;                    
                %imu data accumulator for the covariance update
                imuaccum=imuaccum+[accinc;gyroinc]*dt1;                
            end
            %run strapdown with angle and velocity increments
            [filter.qs2n, filter.Vn, filter.Cen, filter.height]=strapdown_Cen_quat_v000(...
                filter.qs2n, filter.Vn, filter.Cen, filter.height, velinc, angleinc, dt1);     
        end
        function ffun_covariance(filter, imuaccum, covupt_time, curimutime )
            %propagate the covariance corresponds to states, rs in e, vs in e, q s2e,
            % ba, bg, sa, sg
            covdt=curimutime-covupt_time;            
            [STM, Qd]=sys_metric_phipsi_v000(filter.Cen, filter.height, ...
                filter.Vn, filter.qs2n,  imuaccum(1:3)/covdt,imuaccum(4:6)/covdt, ...
                covdt,filter.mode,filter.imuType, filter.imuErrorModel); %use psi implementation with 2
            %              if(filter.imuErrorModel==4)
            %                     [STM, Qd]=sys_metric_phipsi_v002(filter.Cen, filter.height, ...
            %                         filter.Vn, filter.qs2n,  imuaccum(1:3)/covdt,imuaccum(4:6)/covdt, ...
            %                         covdt,filter.mode,filter.imuType); %use psi implementation with 2
            %                     assert(norm(STM1-STM)+norm(Qd1-Qd)<1e-5);
            %              end
            filter.p_k_k=STM*filter.p_k_k*STM'+Qd;  % the covariance of the navigation states and imu error terms
        end
        %==============================================================================================
        %update both the covariance and state
        function correctstates(filter, predict,measure, H,R, ~)
            p_km1_k=filter.p_k_k;
            inno= predict-measure;
            %Kalman
            K=p_km1_k*H'/(H*p_km1_k*H'+R);
            deltaX=K*inno;
            % update covariance
            filter.p_k_k=(eye(filter.covDim)-K*H)*p_km1_k*(eye(filter.covDim)-K*H)'+K*R*K';
            % compute updated state
            [filter.qs2n, filter.Vn, filter.Cen, filter.height]=correctnav_Cen_v000(...
                filter.qs2n, filter.Vn, filter.Cen, filter.height,...
                deltaX(1:filter.imuOrientSIP+2), filter.mode, filter.mechanization);        
            filter.imuErrors = filter.imuErrors + deltaX(filter.imuBiasDriftSIP:end);
        end
        function SaveToFile(filter, inillh_ant, preimutime, ffilres)
            %Write result to the files
            vr_c=rotqr2eu('xyz',filter.qs2n)*180/pi; % Cs2n
            xyz_ant=ecef2geo_v000([asin(-filter.Cen(3,3));asin(-filter.Cen(2,1));filter.height],1)+...
                quatrot_v000(quatmult_v001(rotro2qr(filter.Cen),filter.qs2n,1),filter.Tant2imu,0);
            if(~isempty(inillh_ant))
                vr_a=posdiff_v001(xyz_ant, inillh_ant);
                fwrite(ffilres,[preimutime;vr_a;filter.Vn;vr_c;sqrt(diag(filter.p_k_k(1:9,...
                    1:9)))],'double');
            else
                fwrite(ffilres,[preimutime;ecef2geo_v000(xyz_ant,0);filter.Vn;vr_c;sqrt(diag(filter.p_k_k(1:9,...
                    1:9)))],'double');
            end
        end
        function mag=GetVelocityMag(filter)
            mag=norm(filter.Vn,2); 
        end
        function  SetCamAndRvqs2e(filter, options)
            filter.rvqs2e=zeros(10,1);
            Ce2n=filter.Cen;
            heights2n=filter.height;
            Re=6378137/(sqrt(1.0-0.00669437999014*Ce2n(3,3)^2));
            filter.rvqs2e(1:3,1)=-[(Re+heights2n)*Ce2n(3,1);(Re+heights2n)*Ce2n(3,2);(Re*(1-0.00669437999014)+heights2n)*Ce2n(3,3)];
            filter.rvqs2e(4:6,1)=Ce2n'*filter.Vn;
            filter.rvqs2e(7:10,1)=quatmult_v001(rotro2qr(Ce2n),filter.qs2n,1);            
            filter.camPose=zeros(7,1);
            filter.camPose(5:7)=options.Cimu2cam*options.Cb2imu*(options.Timu2body-options.Tcam2body);
            filter.camPose(1:4)=rotro2qr(options.Cimu2cam);
        end
    end
end