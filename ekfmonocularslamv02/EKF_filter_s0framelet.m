% EKF model for camera imu calibration
% define s0 frame of sensor is tied to epoch t0 of the imu, and c0 frame is
% also tied to this epoch t0, they must be rigidly fixed
% the body frame, in our case is assumed to be aligned with the imu frame
% we assume the s0 frame is well know w.r.t ECEF frame

% the complete state is defined as 
% delta(pos s in s0 frame, vx, vy, vz in s0 frame, q s0 to s frame)
% delta( gravity in s0 frame)
% accelerometer and gyro bias as random walk,
% acc scale and gyro scale factor as random walk, 
% the camera configuration, delta(qimu2cam), and delta(Timu2cam),

% The state dynamics are driven by a 18 (accel and gyro white noise,
% bias white noise and scale factor white noise) dimensional white Gaussian
% noise source and the observations are corrupted by additive white noise
classdef EKF_filter_s0framelet < handle
    properties (Hidden)
        type = 'ekf';
        tag  = 'EKF_CAM_IMU_GRAVITY';  % ID tag
    end
    % The following properties can be set only by class methods
    properties (SetAccess = private)
        invalidateIMUerrors; % reset IMU errors when they goes beyond expectation
        imuType; % the type of IMU, e.g., MEMS 3DX GM 3-35, HG1700, Steval inems, etc.
        imuErrorModel=3;
        % the IMU error model corresponding to random walk bias and
        % scale factor error, 1 without turn on bias estimates,
        % first order GM bias and random walk scale factor errors
        % 2 with random constant turn on bias estimates,
        % first order GM bias and random walk scale factor errors
        % 3 assumes bias and scale factors are random walks
        % 4 assumes random constants bias and scale factor errors
    
        dt; % sampling interval, unit sec    
        rqs02e; % position of s0 in the e frame, qs02e, constant
        rvqs0; % position of current s(k) frame in s0 frame, velocity in s0 frame, q s0 to s(k)    
     
        imuErrors; %ba, bg, sa, sg
        Tant2imu; % antenna position in IMU s-frame
        Cimu2cam;
        camPose; % qs2cam, Ts2cam, i.e., the position of the IMU in the camera frame
      
        imuOrientSIP=7; % imu orientation covariance start index in s-frame
        imuBiasDriftSIP=13; %imu bias error covariance start index
        imuScaleFactorSIP=19;
        camConfigSIP=25;
        covDim=30; % the dimension of P matrix
        p_k_k; % the covariance of the entire state vector 
        camPoseFileHandle;
    end
    methods
        function filter = EKF_filter_s0framelet(options)
            filter.invalidateIMUerrors=options.InvalidateIMUerrors;
            filter.imuType=options.imutype;
            filter.imuErrorModel=options.imuErrorModel;          
            filter.dt=options.dt; 
           
            inillh=options.inillh_ant;
            filter.Tant2imu= options.Cb2imu*(options.Tant2body-options.Timu2body);
            
            initqs2n=quatmult_v001(options.qb2n, rotro2qr(options.Cb2imu), 2); 
            Ce2n0=llh2dcm_v000(inillh(1:2),[0;1]);
            qs2e=quatmult_v001(rotro2qr(Ce2n0), initqs2n,1);           
           
            xyz_imu=ecef2geo_v000(inillh,1)-quatrot_v000(qs2e,filter.Tant2imu,0);            
            filter.rqs02e=[xyz_imu; qs2e];
            filter.rvqs0=zeros(10,1); % start from stationary
            filter.rvqs0(4:6)=options.Vs0;
            filter.rvqs0(7)=1; % identity rotation matrix  
            
            filter.imuErrors=options.imuErrors;
            if(options.useCam)
                filter.camPose=zeros(7,1);
                filter.Cimu2cam=options.Cimu2cam;
                filter.camPose(5:7)=options.Cimu2cam*(-options.Tcam2imu);
                filter.camPose(1:4)=rotro2qr(options.Cimu2cam);
                filter.camPoseFileHandle=fopen(options.camPoseFile,'W');
                % put header
                fprintf(filter.camPoseFileHandle,'%%IMU log time, Cs2c in Euler angles, R,P,Y,X,Y,Z of Ts in c and Cs02e in euler angles, and T s0 in e frame\n');
            else
                filter.camPose=zeros(7,1);
                filter.camPose(1)=1;
            end
            
            filter.p_k_k=zeros(filter.covDim);    % in s0 frame, position, velocity, attitude are of very small variance        
         
            filter.p_k_k(1:3,1:3)=eye(3)*2^2; %position errors in meter
            filter.p_k_k(4:6,4:6)=eye(3)*0.1^2; %vel error in m/s
            filter.p_k_k(filter.imuOrientSIP+(0:2),filter.imuOrientSIP+(0:2))=...
                diag([options.initAttVar,options.initAttVar,options.initAttVar*2].^2);
            filter.p_k_k(10:12,10:12)=diag([1;1;1]*0.02^2); % m/s^2 of gravity in s0 frame
            IMU_ERRDEF=imu_err_defs_v000(options.imutype);
            filter.p_k_k(filter.imuBiasDriftSIP+(0:2),filter.imuBiasDriftSIP+(0:2))=4*eye(3)*IMU_ERRDEF.acc_bias_var;% enlarge initial std by 2
            filter.p_k_k(filter.imuBiasDriftSIP+(3:5),filter.imuBiasDriftSIP+(3:5))=4*eye(3)*IMU_ERRDEF.gyro_bias_var;
            
            filter.p_k_k(filter.imuScaleFactorSIP+(0:2),filter.imuScaleFactorSIP+(0:2))=4*eye(3)*IMU_ERRDEF.acc_scale_var;
            filter.p_k_k(filter.imuScaleFactorSIP+(3:5),filter.imuScaleFactorSIP+(3:5))=4*eye(3)*IMU_ERRDEF.gyro_scale_var;
            filter.p_k_k(25:27,25:27)=diag(([0.5;0.5;0.5]*pi/180).^2); % unit rad
            filter.p_k_k(28:30,28:30)=diag(([4;4;5]*1e-3).^2); % unit m
            %filter.p_k_k=sparse (filter.p_k_k);
        end        
        %===============================================================================================
        %-- State transition function
        % propogate state with accelerometer and gyro input at time k-1 to predict state at k,
        % i.e., X(k|k-1), from state at k-1. U1 contains IMU measurement
        % acc, and gyro angular rate, or delta v, and delta theta
        % and 7th row is the previous epoch(k-1) and 8th is the current
        % epoch(k). Optionally, 9-11th rows are preset gravity in e frame, 
        % 12-14 th row is WIW_E, angular rate of local w-frame w.r.t i-frame
        % represented in e-frame
        % imuaccum records the accumulated delta v and delta angle for two speed
        % covariance update
        function imuaccum= ffun_state(filter, imuaccum, U1, isDelta)
            if(nargin==3) 
                isDelta= false;
            end
            if(isDelta)
                dt1=U1(8,end)-U1(7,1);
                gyroinc=sum(U1(4:6,:),2);
                accinc=sum(U1(1:3,:),2);
                % vr_a=0.5*cross(gyroinc, accinc);
                %accinc=accinc+vr_a; % this line is advised not to use
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
               % vr_a=0.5*cross(gyroinc, accinc)*dt1;
               % accinc=accinc+vr_a; % this line caused worse result on simulation
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
            if(size(U1,1)>8)
                gnomegae=U1(9:end,:);
            else
                %Gravity (most time consuming part of ecef implementations)
                xyz_imu=filter.rqs02e(1:3)+quatrot_v000(filter.rqs02e(4:7),filter.rvqs0(1:3),0);
                Llh=ecef2geo_v000(xyz_imu,0);
                Cn2e=pos2Cne_v000(Llh(1), Llh(2));
                [Rn, Re, gn, sL, cL, WIE_E]=geoparam_v000(Llh);                
                gnomegae=[Cn2e*[0;0;gn];0;0;WIE_E];
            end
            filter.rvqs0=strapdown_local_quat(filter.rvqs0, filter.rqs02e, accinc, gyroinc, dt1, gnomegae);
        end
        function ffun_covariance(filter, imuaccum, covupt_time, curimutime )
            % propagate the covariance in the local s0 frame
            % the covariance corresponds to states, 
            % rs in s0, v s in s0, q s0 2s, gravity in s0, ba, bg, sa, sg, qs2c, Ts in c         
            covdt=curimutime-covupt_time;
            Pvf=filter.p_k_k;
            [STM, Qd]=sys_local_dcm_v001(filter.rqs02e, filter.rvqs0,...
                imuaccum(1:3)/covdt,imuaccum(4:6)/covdt, covdt,filter.imuType, filter.imuErrorModel);
            Pv=Pvf(1:filter.camConfigSIP-1,1:filter.camConfigSIP-1);
            Pv=STM*Pv*STM'+Qd;  % the covariance of the navigation states and imu error terms
            % the total covariance
            Pvf = [ Pv         STM*Pvf(1:filter.camConfigSIP-1,filter.camConfigSIP:end);
                Pvf(filter.camConfigSIP:end,1:filter.camConfigSIP-1)*STM'        Pvf(filter.camConfigSIP:end,filter.camConfigSIP:end)];
            
            filter.p_k_k=Pvf; %sparse(Pvf); % check if Pvf is sparse
        end
        %==============================================================================================
        % if type=0, default, update both the covariance and state
        function correctstates(filter, predict,measure, H,R, type)
            p_km1_k=filter.p_k_k;
            inno= predict-measure;
            %Kalman
            K=p_km1_k*H'/(H*p_km1_k*H'+R);
            deltaX=K*inno;
            % update covariance if desired
            if(nargin==5||type==0)
                filter.p_k_k=(eye(size(p_km1_k,1))-K*H)*p_km1_k*(eye(size(p_km1_k,1))-K*H)'+K*R*K';
                
                % compute updated states                
                filter.rvqs0(1:6)=filter.rvqs0(1:6)-deltaX(1:6); % position and velocity           
                qst=rvec2quat_v000(deltaX(7:9));
                filter.rvqs0(7:10)=quatmult_v000(qst, filter.rvqs0(7:10));                
                filter.imuErrors = filter.imuErrors + deltaX(filter.imuBiasDriftSIP:filter.camConfigSIP-1);
               
                % update camera related states                
                qct=rvec2quat_v000(deltaX(filter.camConfigSIP+(0:2)));
                filter.camPose(1:4) = quatmult_v000(qct, filter.camPose(1:4));
                filter.camPose(5:7) = filter.camPose(5:7)-deltaX(filter.camConfigSIP+(3:5));   
                
            end
        end
        % predict the rotation from c to c0 and translation of c0 in c, and
        % compute the H matrix 
        function [Rc2c0, Tc02c, Hmat]=computeH(filter)
            qrs02c= quatmult_v001( filter.camPose(1:4), filter.rvqs0(7:10) ,0);
            qrc2c0=quatmult_v001(filter.camPose(1:4),qrs02c,2);
            Rc2c0=rotqr2ro(qrc2c0); % predicted Rc2c0 and Tc02c
            Tc02c=filter.camPose(5:7)-quatrot_v000(qrs02c, filter.rvqs0(1:3), 0)-quatrot_v000(qrc2c0, filter.camPose(5:7), 1);
            qrs2c0=quatmult_v001(filter.camPose(1:4), filter.rvqs0(7:10),2);
            
            deltaT=Tc02c-filter.camPose(5:7); 
            % in computing H, use predicted values on the differential point
            % here Rc2c0 and Tc02c are predicted values, RTcnc0 are
            % measured values
            Hmat=sparse([zeros(3,6), -rotqr2ro(qrs2c0), zeros(3,size(filter.p_k_k,1)-15),  eye(3)-Rc2c0, zeros(3,3);
                -rotqr2ro(qrs02c), zeros(3,3), skew(deltaT)*rotqr2ro(filter.camPose(1:4)), ...
                zeros(3,size(filter.p_k_k,1)-15), skew(deltaT)+Rc2c0'*skew(filter.camPose(5:7)), eye(3)-Rc2c0']);
% fixing scale factor and Rs2c
%             Hmat=sparse([zeros(3,6), -rotqr2ro(qrs2c0), zeros(3,size(filter.p_k_k,1)-9);
%                 -rotqr2ro(qrs02c), zeros(3,3), skew(deltaT)*rotqr2ro(filter.camPose(1:4)), ...
%                 zeros(3,size(filter.p_k_k,1)-12), eye(3)-Rc2c0']);
        end      
        % disable the camera configuration states in the state vector
        function disable_camerastates( filter, estimate)
            % store up the estimates
            if(~isempty(estimate))
                filter.camPose=estimate;
            end
            % remove the camera states from the covariance, thus fixing the
            % filter.camPose
            indexFromWichDeleteP=filter.camConfigSIP;
            parToDeleteP=6;
            Pk = [filter.p_k_k(:,1:indexFromWichDeleteP-1) filter.p_k_k(:,indexFromWichDeleteP+parToDeleteP:end)];
            filter.p_k_k =sparse( [Pk(1:indexFromWichDeleteP-1,:); Pk(indexFromWichDeleteP+parToDeleteP:end,:)]);
            
            % set index equal, squeeze camConfig states out in later operations
            filter.groupFrameSIP=filter.camConfigSIP;
        end     
       
        function SaveToFile(filter, inillh_ant, preimutime, ffilres)          
            fwrite(ffilres,[preimutime;filter.rvqs0(1:6);rotqr2eu('xyz', filter.rvqs0(7:10))*180/pi;full(sqrt(diag(filter.p_k_k(1:9,1:9))))],'double');
            
        end
        function SaveCamPoseandRqs02e(filter, imgTime, frmId)
                donkey=unskew(eye(3)-filter.Cimu2cam*(rotqr2ro(filter.camPose(1:4)))')*180/pi; % the deviation from the initial estimate of Cs2c
                mule=rotqr2eu('xyz', filter.rqs02e(4:7))*180/pi; % Cs2e
                fprintf(filter.camPoseFileHandle,...
                    '%d\t%15.8f\t%10.6f\t%10.6f\t%10.6f\t%8.5f\t%8.5f\t%8.5f\t%10.6f\t%10.6f\t%10.6f\t%8.5f\t%8.5f\t%8.5f\n',...
                    frmId, imgTime, donkey, quatrot_v000(filter.camPose(1:4),-filter.camPose(5:7),1),mule, filter.rqs02e(1:3));           
        end
     
        function mag=GetVelocityMag(filter)
            mag=norm(filter.rvqs0(4:6),2);
        end      
        function qs2e= get_qs2e(filter)
            qs2e= quatmult_v001(filter.rqs02e(4:7), filter.rvqs0(7:10),2);
        end
        function rsine= get_rs2e(filter)
            rsine= filter.rqs02e(1:3) +quatrot_v000(filter.rqs02e(4:7),filter.rvqs0(1:3),0);
        end
    end
end