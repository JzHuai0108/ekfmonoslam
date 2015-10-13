% use EKF to estimate gyro and accel biases, scale factor is ignored, for
% MEMS IMU grade close to microstrain 3dm gx3-35
% define s0 frame of sensor is tied to w-frame(n-frame rotated about z-axis) at epoch t0

% the body frame, in our case is assumed to be aligned with the imu frame
% assume constant gravity in s0 frame, which is estimated as the average of
% accelerometer measurements in the static session
% the complete state is defined as 
% delta(pos s in s0 frame, vx, vy, vz in s0 frame, q s0 to s frame)
% accelerometer and gyro bias as random walk, or constants

% attitude expressed in quaternion, taking 4 dimensions, but in the P
% matrix, their covariance are in terms of radians^2, taking only 3 dimensions.
% The P matrix is in sparse form, not full form.

classdef EKF_filter_s0frame_bias < handle
    properties (Hidden)
        type = 'ekf';
        tag  = 'EKF_GPS_IMU_GRAVITY';  % ID tag
    end
    % The following properties can be set only by class methods
    properties (SetAccess = private)
        invalidateIMUerrors; % reset IMU errors when they goes beyond expectation
        imuType; % the type of IMU, e.g., MEMS 3DX GM 3-35, Steval inems, etc.
        imuErrorModel=3;
        % the IMU error model corresponding to bias
        dt; % sampling interval, unit sec        
        rqs02e; % position of s0 in the e frame, qs02e, constant
        rvqs0; % position of current s(k) frame in s0 frame, velocity in s0 frame, q s0 to s(k)    
        imuErrors; %ba, bg     
        imuOrientSIP=7; % imu orientation covariance start index in s-frame
        imuBiasDriftSIP=10; %imu bias error covariance start index      
        covDim=15; % the dimension of P matrix
        velNoiseStd; % noise density driving horizontal velocity
        p_k_k; % the covariance of the entire state vector
    end
    methods
        function filter = EKF_filter_s0frame_bias(options)
            filter.invalidateIMUerrors=options.InvalidateIMUerrors;
            filter.imuType=options.imutype;
            filter.imuErrorModel=options.imuErrorModel;          
            filter.dt=options.dt; 
            
            inillh=options.inillh_ant;            
            Ce2n0=llh2dcm_v000(inillh(1:2),[0;1]);           
            filter.rqs02e=[ecef2geo_v000(options.inillh_ant,1); rotro2qr(Ce2n0')];
            % assume the IMU is very close to antenna
            filter.rvqs0=zeros(10,1);
            filter.rvqs0(7:10)=[options.qb2n(1);-options.qb2n(2:4)]; 
            filter.rvqs0(4:6)=options.Vs0;
           
            filter.imuErrors=options.imuErrors;   
            filter.velNoiseStd=options.velNoiseStd;       
            filter.p_k_k=zeros(filter.covDim);    % in s0 frame, position, velocity, attitude are of very small variance        
            
            filter.p_k_k(1:3,1:3)=diag([0.2, 0.2, 0.5].^2); %position errors in meter
            filter.p_k_k(4:6,4:6)=diag([0.1, 0.1, 0.5].^2); %vel error in m/s
            filter.p_k_k(filter.imuOrientSIP+(0:2),filter.imuOrientSIP+(0:2))=...
                diag(([options.initAttVar,options.initAttVar,options.initAttVar*3]).^2);
            IMU_ERRDEF=imu_err_defs_v000(options.imutype);
            filter.p_k_k(filter.imuBiasDriftSIP+(0:2),filter.imuBiasDriftSIP+(0:2))=4*eye(3)*IMU_ERRDEF.acc_bias_var;% enlarge initial std by 2
            filter.p_k_k(filter.imuBiasDriftSIP+(3:5),filter.imuBiasDriftSIP+(3:5))=4*eye(3)*IMU_ERRDEF.gyro_bias_var;    
        
            filter.p_k_k=sparse (filter.p_k_k);
        end        
        %===============================================================================================
        %-- State transition function
        % propogate state with accelerometer and gyro input at time k-1 to predict state at k,
        % i.e., X(k|k-1), from state at k-1. U1 contains IMU measurement
        % acc, and gyro angular rate, or delta v, and delta theta
        % and 7th row is the previous epoch(k-1) and 8th is the current
        % epoch(k). Optionally, 9-11th rows are preset gravity in s0 frame, 
        % 12-14 th row is $\omega_{ie}_{s_0}$, angular rate of local w-frame w.r.t i-frame
        % represented in w-frame
        % imuaccum records the accumulated delta v and delta angle for two speed
        % covariance update
        function imuaccum= ffun_state(filter, imuaccum, U1, isDelta, isConstantVel)
            if(isDelta)
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
                        filter.imuErrors=zeros(6,1);
                    end
                end                
                angleinc=gyroinc-filter.imuErrors(4:6)*dt1;
                velinc=accinc-filter.imuErrors(1:3)*dt1;
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
                        filter.imuErrors=zeros(6,1);
                    end
                end
                gyroinc=gyroinc-filter.imuErrors(4:6);
                accinc=accinc-filter.imuErrors(1:3);
                %imu data accumulator for the covariance update
                imuaccum=imuaccum+[accinc;gyroinc]*dt1;
            end           
            if(size(U1,1)>8)
                gwomegaw=U1(9:end,:);
            else
                %Gravity (most time consuming part of ecef implementations)
                xyz_imu=filter.rqs02e(1:3)+quatrot_v000(filter.rqs02e(4:7),filter.rvqs0(1:3),0);
                Llh=ecef2geo_v000(xyz_imu,0);
                Cn2e=pos2Cne_v000(Llh(1), Llh(2));
                [Rn, Re, gn, sL, cL, WIE_E]=geoparam_v000(Llh);                
                gwomegaw=[quatrot_v000(filter.rqs02e(4:7), Cn2e*[0;0;gn], 1);quatrot_v000(filter.rqs02e(4:7),[0;0;WIE_E],1)];
            end
            if(~isConstantVel)
                filter.rvqs0=strapdown_local_quat_bias(filter.rvqs0, filter.rqs02e, accinc, gyroinc, dt1, gwomegaw);    
            else
                filter.rvqs0=strapdown_local_quat_constvel(filter.rvqs0, filter.rqs02e, accinc, gyroinc, dt1, gwomegaw);   
            end
        end
        function ffun_covariance(filter, imuaccum, covupt_time, curimutime, isConstantVel )
            % propagate the covariance in the local s0 frame
            % the covariance corresponds to states, 
            % rs in s0, v s in s0, q s0 2s, gravity in s0, ba, bg, sa, sg, qs2c, Ts in c,
            % for each group frame, qs02si and Tsi in s0, for each
            % point, its inverse depth
            covdt=curimutime-covupt_time;     
            if(~isConstantVel)
            [STM, Qd]=sys_local_dcm_bias(filter.rqs02e, filter.rvqs0,...
                imuaccum(1:3)/covdt,imuaccum(4:6)/covdt, covdt,filter.imuType, filter.imuErrorModel); 
            else
                [STM, Qd]=sys_local_dcm_constvel(filter.rqs02e, filter.rvqs0,...
                imuaccum(1:3)/covdt,imuaccum(4:6)/covdt, covdt,filter.imuType, filter.velNoiseStd); 
            end
            % the total covariance
            Pvf = STM*filter.p_k_k*STM'+Qd;  % the covariance of the navigation states and imu error terms
            filter.p_k_k=sparse(Pvf); % check if Pvf is sparse
        end
        %==============================================================================================
        % if type==0, default, update all the state
        % if type==1, update only position, velocity, accelerometer bias
        % this trick is inspired by S. Weiss PhD diss. page 77, Sec 3.1.6
        function correctstates(filter, predict,measure, H,R, type)
            p_km1_k=filter.p_k_k;
            inno= predict-measure;
            %Kalman gain
            K=p_km1_k*H'/(H*p_km1_k*H'+R);
            deltaX=K*inno;            
            if(nargin==5||type==0)
                filter.p_k_k=(eye(size(p_km1_k,1))-K*H)*p_km1_k*(eye(size(p_km1_k,1))-K*H)'+K*R*K';
                % compute updated states                
                filter.rvqs0(1:6)=filter.rvqs0(1:6)-deltaX(1:6); % position and velocity           
                qst=rvec2quat_v000(deltaX(7:9));
                filter.rvqs0(7:10)=quatmult_v000(qst, filter.rvqs0(7:10));                
                filter.imuErrors = filter.imuErrors + deltaX(filter.imuBiasDriftSIP:end);    
            elseif(type==1)
                filter.p_k_k=(eye(size(p_km1_k,1))-K*H)*p_km1_k*(eye(size(p_km1_k,1))-K*H)'+K*R*K';
                % update states                
                filter.rvqs0(1:6)=filter.rvqs0(1:6)-deltaX(1:6); % position and velocity                            
                filter.imuErrors(1:3) = filter.imuErrors(1:3) + deltaX(filter.imuBiasDriftSIP+(0:2));  
            end
        end      
       
        function SaveToFile(filter, inillh_ant, preimutime, ffilres)          
            fwrite(ffilres,[preimutime;filter.rvqs0(1:6);rotqr2eu('xyz', [filter.rvqs0(7); -filter.rvqs0(8:10)])*180/pi;full(sqrt(diag(filter.p_k_k(1:9,1:9))))],'double');
        end  
        function mag=GetVelocityMag(filter)
            mag=norm(filter.rvqs0(4:6),2);
        end  

    end
end