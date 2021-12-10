% EKF model for camera gps imu integration
% The following state space model is used :
% the complete state is defined as delta(ecef x,y,z of sensor at epoch t0,
% q s0 2 e frame(s0 frame of sensor is tied to epoch t0 of the imu)),
% the delta(pos s in s0 frame, vx, vy, vz in s0 frame of sensor, q s0 to s frame)
% accelerometer and gyro bias as random walk,
% acc scale and gyro scale factor as random walk,
% the camera configuration, delta(qimu2cam), and delta(Timu2cam),
% follows by delta of group frames and feature invese depths.

% the camera configuration can be turned off after they
% ramain constant for a while, and their values are stored in model.
% attitude expressed in quaternion, taking 4 dimensions, but in the P
% matrix, their covariance are in terms of radians^2, taking only 3 dimensions.
% The P matrix is in sparse form, not full form.
% The state dynamics are driven by a 18 (accel and gyro white noise,
% bias white noise and scale factor white noise) dimensional white Gaussian
% noise source and the observations are corrupted by additive white noise
classdef EKF_filter_s0frame < handle
    properties (Hidden)
        type = 'ekf';
        tag  = 'EKF_CAM_IMU_GPS';  % ID tag
    end
    % The following properties can be set only by class methods
    properties (SetAccess = private)
        invalidateIMUerrors; % reset IMU errors when they goes beyond expectation
        imuType; % the type of IMU, e.g., MEMS 3DX GM 3-35, HG1700, etc.
        imuErrorModel=3;
        % the IMU error model corresponding to random walk bias and
        % scale factor error, 1 without turn on bias estimates,
        % first order GM bias and random walk scale factor errors
        % 2 with random constant turn on bias estimates,
        % first order GM bias and random walk scale factor errors
        % 3 assumes bias and scale factors are random walks
        % 4 assumes random constants bias and scale factor errors
        % the body frame, in our case is assumed to be aligned with the H764G frame,
        % but here we are testing on an arbitray sensor, like a mems imu, so Cb2imu
        % represent the rotations between the body frame and the imu frame
        % set it as [] to be determined using IMU coarse alignment
        Cb2imu;
        dt; % sampling interval, unit sec
        % the translation from antenna to mems frame, i.e., the antenna's
        % position in the mems frame.
        Tant2imu;
        
        camType;
        sigmaCAM; % camera measurement std dev
        rvqs2e; % position of IMU in e frame, v in e, and qs2e
        rqs02e; % position of s0 in the e frame, qs02e
        rvqs0; % position of s frame in s0 frame, velocity in s0 frame, q s0 to s
        imuErrors; %ba, bg, sa, sg
        camPose; % qs2cam, Ts2cam, i.e., the position of the IMU in the camera frame
        groupPose; % the poses of all the group frames, it is a group mapping list
        % all the poses are in the states
        % the qs0 2 si, and Tsi in s0, the attitude and position of group frame si in the s0 frame
        % each element of the groupPose list, has a unique group number,
        % number of points, and qs02si and Tsi in s0
        features_info; % the feature points list in the states vector
        % each element of this list has many attributes about the feature point
        imuqs02eSIP=4; % imu orientation covariance start index in e-frame
        imuOrientSIP=13; % imu orientation covariance start index in s-frame
        imuBiasDriftSIP=16; %imu bias error covariance start index
        imuScaleFactorSIP=22;
        camConfigSIP=28;
        groupFrameSIP=34; % the covariance start index of group frame orientation and positions
        featSIP=34;
        refGrpFrmId=1;
        p_k_k; % the covariance of the entire state vector
        initial_rho=0.1;
        std_rho=0.1;
        camPoseFileHandle;
    end
    methods
        % featPts is a feature list, each column is a feature point structure
        % heritage is the previous filter, contains all details of
        % GPS IMU navigation states
        function filter = EKF_filter_s0frame(options, heritage, featPts)
            filter.invalidateIMUerrors=options.InvalidateIMUerrors;
            filter.imuType=options.imutype;
            filter.imuErrorModel=heritage.imuErrorModel;
            filter.camType=options.camtype;
            filter.dt=options.dt;
            filter.sigmaCAM=options.sigmaCAM;
            filter.Cb2imu=options.Cb2imu;      
            filter.Tant2imu=filter.Cb2imu*(options.Tant2body-options.Timu2body);        
      
            % these two values in camPose, Timu2cam and Cimu2cam of the 
            % camera configuration are to be refined in the filtering process and then
            % frozen after they remain almost constant with enough image
            filter.camPose=heritage.camPose;     
            filter.camPoseFileHandle=fopen(options.camPoseFile,'W');
            % put header
            fprintf(filter.camPoseFileHandle,'%%IMU log time, Cs2c in Euler angles, R,P,Y,X,Y,Z of Ts in c\n');
            % Initial state vector and covariance matrix
            filter.rvqs2e= heritage.rvqs2e;
            filter.rqs02e=heritage.rvqs2e([1:3, 7:10]);
            filter.rvqs0=zeros(10,1);
            filter.rvqs0(4:6)=quatrot_v000(heritage.rvqs2e(7:10),heritage.rvqs2e(4:6),1);
            filter.rvqs0(7)=1; % identity rotation matrix
            filter.imuErrors=heritage.imuErrors;
            
            grpPose.grpId=featPts(1).grpId;
            grpPose.grpFrmNo=featPts(1).initFrmNo;
            grpPose.pose=filter.rvqs0([7:10,1:3]); % for the first group frame
            grpPose.ptsNo=length(featPts);
            % for groupPose and features_info column corresponds to observations and row corresponds to attributes
            filter.groupPose=grpPose;
            filter.features_info=featPts;
            
            covDim=filter.groupFrameSIP-1+length(filter.groupPose)*6+length(filter.features_info);
            filter.p_k_k=zeros(covDim);
            Gmat=zeros((2+3+4)*3, (3+4)*3);
            Gmat(1:3,1:3)=eye(3);
            Gmat(4:6,7:9)=eye(3);
            Cs02e=quat2dcm_v000(filter.rvqs2e(7:10));
            Gmat(10:12,4:9)=[Cs02e', Cs02e'*skew(-filter.rvqs2e(4:6))];
            Gmat(16:end,10:end)=eye(12);
            filter.p_k_k(1:(2+3+4)*3, 1:(2+3+4)*3)=Gmat*heritage.p_k_k*Gmat';
            
            filter.p_k_k(filter.camConfigSIP+(0:2),filter.camConfigSIP+(0:2))=diag([options.initAttVar,options.initAttVar,options.initAttVar*2].^2);
            filter.p_k_k(filter.camConfigSIP+(3:5),filter.camConfigSIP+(3:5))=diag([0.1;0.1;0.1].^2);
            filter.refGrpFrmId=grpPose.grpId;
            % the reference group frame has very small variance
            smallValue=0;
%             assert(sum(sum(filter.p_k_k(filter.imuOrientSIP+(0:2), filter.imuOrientSIP+(0:2))))+...
%                 sum(sum(filter.p_k_k(7:9, 7:9)))==0);
            filter.p_k_k(filter.groupFrameSIP+(0:2),filter.groupFrameSIP+(0:2))=smallValue;
            filter.p_k_k(filter.groupFrameSIP+(3:5),filter.groupFrameSIP+(3:5))=smallValue;
            
            filter.featSIP=filter.groupFrameSIP+length(filter.groupPose)*6;
            filter.p_k_k(filter.featSIP:end, filter.featSIP:end)=eye(length(filter.features_info))*filter.std_rho^2;
            filter.p_k_k=sparse (filter.p_k_k);
        end
        
        %===============================================================================================
        %-- State transition function
        % propogate state with accelerometer and gyro input at time k-1 to predict state at k,
        % i.e., X(k|k-1), from state at k-1. U1 contains IMU measurement
        % acc, and gyro angular rate
        % and 7th row is the previous epoch(k-1) and 8th is the current epoch(k)
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
            
            filter.rvqs0=[quatrot_v000(filter.rqs02e(4:7),filter.rvqs2e(1:3)-filter.rqs02e(1:3),1);...
                quatrot_v000(filter.rqs02e(4:7),filter.rvqs2e(4:6),1);...
                quatmult_v001(filter.rvqs2e(7:10),filter.rqs02e(4:7),1)];
            
        end
        function ffun_covariance(filter, imuaccum, covupt_time, curimutime)
            %propagate the covariance in the local s0 frame
            % the covariance corresponds to states, rs0 in e, q s02e,
            % rs in s0, v s in s0, q s02s, ba, bg, sa, sg, qs2c, Ts in c,
            % for each group frame, qs02si and Tsi in s0, for each
            % point, its inverse depth
            covdt=curimutime-covupt_time;
            Pvf=filter.p_k_k;
            [STM, Qd]=sys_local_dcm_v000(filter.rqs02e, filter.rvqs0,...
                imuaccum(1:3)/covdt,imuaccum(4:6)/covdt, covdt,filter.imuType, filter.imuErrorModel);
            Pv=Pvf(1:filter.camConfigSIP-1,1:filter.camConfigSIP-1);
            Pv=STM*Pv*STM'+Qd;  % the covariance of the navigation states and imu error terms
            % the total covariance
            Pvf = [ Pv         STM*Pvf(1:filter.camConfigSIP-1,filter.camConfigSIP:end);
                Pvf(filter.camConfigSIP:end,1:filter.camConfigSIP-1)*STM'        Pvf(filter.camConfigSIP:end,filter.camConfigSIP:end)];
            
            % Propagate the imu error states, the camera configuration, and the
            % feature states, in general make no difference so usually commented
            assert(sum(sum(STM(filter.imuBiasDriftSIP:filter.camConfigSIP-1,filter.imuBiasDriftSIP:filter.camConfigSIP-1)))==12);
            %             filter.imuErrors=STM(filter.imuBiasDriftSIP:filter.camConfigSIP-1,filter.imuBiasDriftSIP:filter.camConfigSIP-1)*filter.imuErrors;
            filter.p_k_k=Pvf;
        end
        %==============================================================================================
        % if type=0, default, update both the covariance and state, if type=1,
        % only return corrected states, X_k_k, the filter states remain intact,
        % for 1 point ransac, X_k_k is all the states connected
        function X_k_k=correctstates(filter, predict,measure, H,R, type)
            p_km1_k=filter.p_k_k;
            inno= predict-measure;
            %Kalman
            K=p_km1_k*H'/(H*p_km1_k*H'+R);
            deltaX=K*inno;
            % update covariance if desired
            if(nargin==5||type==0)
                filter.p_k_k=(eye(size(p_km1_k,1))-K*H)*p_km1_k*(eye(size(p_km1_k,1))-K*H)'+K*R*K';
            end
            % compute updated state
            if(nargin<6||type==0)
                [filter.rqs02e, filter.rvqs0]=...
                    correctnav_localframe_v000(filter.rqs02e, filter.rvqs0, deltaX(1:filter.imuOrientSIP+2));
                % reflect the update in the global frame
                filter.rvqs2e(1:3)=filter.rqs02e(1:3)+quatrot_v000(filter.rqs02e(4:7),filter.rvqs0(1:3),0);
                filter.rvqs2e(4:6)=quatrot_v000(filter.rqs02e(4:7),filter.rvqs0(4:6),0);
                filter.rvqs2e(7:10)=quatmult_v001(filter.rqs02e(4:7),filter.rvqs0(7:10),2);
                
                filter.imuErrors = filter.imuErrors + deltaX(filter.imuBiasDriftSIP:filter.camConfigSIP-1);
                % update camera related states
                if(filter.camConfigSIP~=filter.groupFrameSIP)
                    % camera configuration states still in the states
                    qct=rvec2quat_v000(deltaX(filter.camConfigSIP+(0:2)));
                    filter.camPose(1:4) = quatmult_v000(qct, filter.camPose(1:4));
                    filter.camPose(5:7) = filter.camPose(5:7)-deltaX(filter.camConfigSIP+(3:5));
                end
                for snag=1: length(filter.groupPose)
                    qsi2t=rvec2quat_v000(deltaX(filter.groupFrameSIP+(snag-1)*6+(0:2))); % si frame to true si frame
                    filter.groupPose(snag).pose(1:4) = quatmult_v000(qsi2t, filter.groupPose(snag).pose(1:4));
                    filter.groupPose(snag).pose(5:7) = filter.groupPose(snag).pose(5:7)-deltaX(filter.groupFrameSIP+(snag-1)*6+(3:5));
                end
                filter.featSIP=filter.groupFrameSIP+length(filter.groupPose)*6;
                for smug=1: length(filter.features_info)
                    filter.features_info(smug).invDepth=filter.features_info(smug).invDepth-...
                        deltaX(filter.featSIP+smug-1);
                end
            else % implemented for 1-Point RANSAC
                if(filter.camConfigSIP~=filter.groupFrameSIP)
                    % 3 for the 3 4 vectors of quaternions
                    navDim=3*11+3;
                else
                    navDim=3*9+2;
                end
                stateDim=navDim+length(filter.groupPose)*7+length(filter.features_info);
                X_k_k=zeros(stateDim,1);
                [X_k_k(1:7), X_k_k(8:17)]=...
                    correctnav_localframe_v000(filter.rqs02e, filter.rvqs0, deltaX(1:filter.imuOrientSIP+2));
                
                X_k_k(17+(1:12)) = filter.imuErrors + deltaX(filter.imuBiasDriftSIP:filter.camConfigSIP-1);
                % update camera related states
                if(filter.camConfigSIP~=filter.groupFrameSIP)
                    % camera configuration states still in the states
                    qct=rvec2quat_v000(deltaX(filter.camConfigSIP+(0:2)));
                    X_k_k(navDim+(-6:-3)) = quatmult_v000(qct, filter.camPose(1:4));
                    X_k_k(navDim+(-2:0)) = filter.camPose(5:7)-deltaX(filter.camConfigSIP+(3:5));
                end
                
                for snag=1: length(filter.groupPose)
                    qsi2t=rvec2quat_v000(deltaX(filter.groupFrameSIP+(snag-1)*6+(0:2))); % si frame to true si frame
                    X_k_k(navDim+(snag-1)*7+(1:4)) = quatmult_v000(qsi2t, filter.groupPose(snag).pose(1:4));
                    X_k_k(navDim+(snag-1)*7+(5:7)) = filter.groupPose(snag).pose(5:7)-deltaX(filter.groupFrameSIP+(snag-1)*6+(3:5));
                end
                groupDim=navDim+snag*7;
                filter.featSIP=filter.groupFrameSIP+length(filter.groupPose)*6;
                for smug=1: length(filter.features_info)
                    X_k_k(groupDim+smug)=filter.features_info(smug).invDepth-...
                        deltaX(filter.featSIP+smug-1);
                end
            end
        end
        % predict the image cooridnates of points in the states, and
        % compute their H matrix
        function predictPoints(filter, cam)
            % for each group frame cj and current camera frame c, compute Ccj2c and Tcj2c
            grpLeg=length(filter.groupPose);
            groupTransform=zeros(7, grpLeg);
            qs02c=quatmult_v001(filter.camPose(1:4),filter.rvqs0(7:10),0);
            for apple=1:grpLeg
                qs02cj=quatmult_v001(filter.camPose(1:4),filter.groupPose(apple).pose(1:4),0);
                groupTransform(1:4,apple)=quatmult_v001(qs02c,qs02cj,2);
                groupTransform(5:7,apple)=-quatrot_v000(groupTransform(1:4,apple), filter.camPose(5:7),0)+...
                    quatrot_v000(qs02c,filter.groupPose(apple).pose(5:7)-filter.rvqs0(1:3),0)+...
                    filter.camPose(5:7);
            end
            groupIds=[filter.groupPose.grpId];
            % the dimension of the P covariance matrix
            covDim=filter.groupFrameSIP-1+length(filter.groupPose)*6+length(filter.features_info);
            %for each feature with lostNo==0, not lost, predict its position
            for wasp = 1:length(filter.features_info)
                if (filter.features_info(wasp).lostNo==0)
                    % the feature is still being tracked, note some point
                    % may be lost. But before a new group comes, it is still in
                    % the states
                    apple=find(groupIds==filter.features_info(wasp).grpId);
                    % grpId may contain gap, this instruction can be sped up by bisection search
                    rho=filter.features_info(wasp).invDepth;
                    hc=quatrot_v000(groupTransform(1:4,apple), [filter.features_info(wasp).xyn; 1],0)+...
                        rho*groupTransform(5:7,apple); % homogeneous coordiantes in the current camera frame
                    runner=hi_inverse_depth_v002(hc, cam); % project onto the image plane
                    filter.features_info(wasp).h=runner;
                    if ~isempty(runner)
                        % compute the H matrix for this candidate point, this
                        % operation can be optimized a little by recasting the
                        % code, to compute H only for points having measurement
                        dhrl_dx = zeros(3, covDim);
                        hcx = hc(1);
                        hcy = hc(2);
                        hcz = hc(3);
                        xyn=hc(1:2)/hc(3);
                        dz_dxydistort=[cam.fc(1), cam.alpha*cam.fc(1);0, cam.fc(2)];
                        dxyd_dxyn=xydistort_by_xynormalized(xyn, cam);
                        dxyn_dhrl=[1/(hcz)       0           -hcx/(hcz^2);
                            0               1/(hcz)    -hcy/(hcz^2)];
                        Cs02c=quat2dcm_v000(qs02c);
                        dhrl_dx(:,7:9)=-rho*Cs02c;
                        Cs2c=quat2dcm_v000(filter.camPose(1:4));
                        dhrl_dx(:,filter.imuOrientSIP+(0:2))=skew(hc-rho*filter.camPose(5:7))*Cs2c;
                        Ccj2c=quat2dcm_v000(groupTransform(1:4,apple));
                        A=[filter.features_info(wasp).xyn; 1]-rho*filter.camPose(5:7);
                        % derivative w.r.t inverse depth
                        index_of_insertion =filter.groupFrameSIP-1+length(filter.groupPose)*6+wasp;
                        dhrl_dx(:,index_of_insertion)=groupTransform(5:7,apple);
%                         if(filter.features_info(wasp).grpId~=filter.refGrpFrmId)
% u may disable this if line if u set the reference frame variance to very small values
                            % the feature point is not in the reference group frame
                            index_of_insertion =filter.groupFrameSIP+(apple-1)*6;
                            dhrl_dx(:,index_of_insertion+(0:5))=[-Ccj2c*skew(A)*Cs2c, rho*Cs02c];
%                         end
                        % still need to estimate camera configurations
                        if(filter.camConfigSIP~=filter.groupFrameSIP)
                            dhrl_dx(:,filter.camConfigSIP+(0:5))=[skew(hc-rho*filter.camPose(5:7))-Ccj2c*skew(A),rho*(eye(3)-Ccj2c)];
                        end
                        filter.features_info(wasp).H =sparse(dz_dxydistort*dxyd_dxyn*dxyn_dhrl*dhrl_dx);
                        filter.features_info(wasp).S = filter.features_info(wasp).H*filter.p_k_k*filter.features_info(wasp).H' + eye(2)*filter.sigmaCAM^2;
                    else
                        filter.features_info(wasp).lostNo=1;
                    end
                else
                    filter.features_info(wasp).lostNo=filter.features_info(wasp).lostNo+1;
                    filter.features_info(wasp).h =[];
                    filter.features_info(wasp).H=[];
                    filter.features_info(wasp).S=[];
                end
            end
        end
        % determine the number of support points for each hypothesis
        % xi, the hypothesis states computed given only one point's measurement to KF filter
        % cam, the camera intrinsics and distortions
        % z_id, all the measurements of features expressed in inverse depth
        % features_info, the feature points in the states,
        % output: hypothesis_support, number of support points,
        % positions_li_inliers_id, the index of inliers w.r.t all points having
        % observations, i.e., z~=[], this why we need map it
        % to the indices of all points
        function [hypothesis_support, positions_li_inliers_id] = compute_hypothesis_support_fast(filter, xi, cam, threshold )
            hypothesis_support = 0;
            z_id=[filter.features_info.z]; % get all the observations, 2 x N matrix
            h_distorted=zeros(size(z_id));
            % debug
            probe=sum([filter.features_info.lostNo]==0);
            if(probe~=size(z_id,2))
                fprintf('Error in compute_hypothesis_support_fast_v002, features_info.lostNo does not corresponds to z!\n');
            end
            if ~isempty(z_id)
                % for each group frame, compute updated Ccj2c and Tcj2c
                % based on xi
                grpLeg=length(filter.groupPose);
                groupTransform=zeros(7, length(filter.groupPose));
                
                if(filter.camConfigSIP~=filter.groupFrameSIP)
                    % 3 for the 3 4 vectors of quaternions
                    navDim=3*11+3;
                    qs02c=quatmult_v001(xi(navDim+(-6:-3)),xi(14:17),0);
                else
                    navDim=3*9+2;
                    qs02c=quatmult_v001(filter.camPose(1:4),xi(14:17),0);
                end
                groupDim=navDim+length(filter.groupPose)*7;
                if(filter.camConfigSIP~=filter.groupFrameSIP)
                    for apple=1:grpLeg
                        qs02cj=quatmult_v001(xi(navDim+(-6:-3)),xi(navDim+(apple-1)*7+(1:4)),0);
                        groupTransform(1:4,apple)=quatmult_v001(qs02c,qs02cj,2);
                        groupTransform(5:7,apple)=-quatrot_v000(groupTransform(1:4,apple), xi(navDim+(-2:0)),0)+...
                            quatrot_v000(qs02c,xi(navDim+(apple-1)*7+(5:7))-xi(8:10),0)+...
                            xi(navDim+(-2:0));
                    end
                else
                    for apple=1:grpLeg
                        qs02cj=quatmult_v001(filter.camPose(1:4),xi(navDim+(apple-1)*7+(1:4)),0);
                        groupTransform(1:4,apple)=quatmult_v001(qs02c,qs02cj,2);
                        groupTransform(5:7,apple)=-quatrot_v000(groupTransform(1:4, apple), filter.camPose(5:7),0)+...
                            quatrot_v000(qs02c,xi(navDim+(apple-1)*7+(5:7))-xi(8:10),0)+...
                            filter.camPose(5:7);
                    end
                end
                %for each feature with ~isempty(z), as this
                % function is called after the tracker that updates the
                % features_info, predict its position with current hypotheses
                zIter=0; % count how many points have measurements
                groupIds=[filter.groupPose.grpId];
                for wasp = 1:length(filter.features_info)
                    if (~isempty(filter.features_info(wasp).z)) % this is equivalent to features_info(wasp).lostNo==0
                        % the feature is still being tracked, note some point
                        % may be lost. But before a new group comes, it is still in
                        % the states
                        apple=find(groupIds==filter.features_info(wasp).grpId,1);
                        % because grpId contains gap and monotonic increasing,
                        % this instruction can be sped up by making use of bisection search
                        rho=xi(groupDim+wasp);
                        hc=quatrot_v000(groupTransform(1:4,apple), [filter.features_info(wasp).xyn; 1],0)+...
                            rho*groupTransform(5:7,apple);
                        zIter=zIter+1;
                        % we cannot use hi_inverse_depth_v002 because there
                        % may be points with z near the border but with illegitimate uv_d
                        uv_d=project_points2(hc,zeros(3,1),zeros(3,1),cam.fc,cam.cc,cam.kc,cam.alpha); %uv_d 2xN point vector in pixels
                        h_distorted(:,zIter)=uv_d;                   
                    end
                end
                nu = z_id - h_distorted;
                residuals = sqrt(nu(1,:).^2+nu(2,:).^2);
                positions_li_inliers_id = residuals<threshold;
                hypothesis_support = hypothesis_support + sum(positions_li_inliers_id);
            else
                positions_li_inliers_id = [];
            end
        end
        
        function ransac_hypotheses( filter, cam )
            p_at_least_one_spurious_free = 0.99; % default value
            % RANSAC threshold should have a low value (less than the standard
            % deviation of the filter measurement noise); as high innovation points
            % will be later rescued
            % Klein and Murray also use 1 for 0 pyramid level noise
            % and my preliminary tests show changing from 1 sigma to 2
            % sigma does not improve performance
            threshold = filter.sigmaCAM;             
            max_hypothesis_support = 0; % will be updated
            individually_compatible =[filter.features_info.lostNo]==0; % indicate if a point has been tracked        
            % debug
            for strut=1:length(filter.features_info)
                if(filter.features_info(strut).lostNo==0)
                    assert(~isempty(filter.features_info(strut).z));
                end
            end
            % debug end
            num_IC_matches = sum(individually_compatible);
            n_hyp = min(num_IC_matches*2, num_IC_matches+5); % initial number of iterations, will
            % be updated if at least one point is an inlier
            % generate a mapping table that maps the index of zi
            % (each measurement) to the index of features_info
            zmap=find(individually_compatible);
            best_li_inliers=[];
            for i = 1:n_hyp
                % select random match
                random_match_position = floor(rand(1)*num_IC_matches)+1;
                position = zmap(random_match_position);
                zi = filter.features_info(position).z;
                % 1-match EKF state update
                hi = filter.features_info(position).h;             
                Hi = sparse(filter.features_info(position).H);
                % do not update the states and covariance of the filter,
                % only return the corrected states as Civera, 2010.
                xi= filter.correctstates(hi,zi, Hi, eye(2)*filter.sigmaCAM^2, 1);
                
                % Compute hypothesis support: predict measurements and count matches
                % under a threshold
                [hypothesis_support, positions_li_inliers_id] ...
                    = filter.compute_hypothesis_support_fast(xi, cam, threshold );
                
                if hypothesis_support > max_hypothesis_support
                    max_hypothesis_support = hypothesis_support;
                    best_li_inliers=positions_li_inliers_id;
                    epsilon = 1-(hypothesis_support/num_IC_matches);
                    n_hyp = ceil(log(1-p_at_least_one_spurious_free)/log(1-(1-epsilon)));
                    if n_hyp==0
                        break;
                    end
                end
                if i>n_hyp
                    break;
                end
            end
            ftPtId=zmap(best_li_inliers);
            for smile=1:length(ftPtId)
                filter.features_info(ftPtId(smile)).low_innovation_inlier = 1;
            end
        end
        % update KF based on low innovation inliers
        function ekf_update_li_inliers( filter )
            % mount vectors and matrices for the update
            allocsz=sum([filter.features_info.low_innovation_inlier]);
            if(allocsz==0)
                fprintf('No points in the states are low inliers!\n');
                return;
            end
            z=zeros(allocsz*2,1);
            h=z;
            H=zeros(allocsz*2, size(filter.p_k_k,1));
            featurectr=0;
            for i=1:length( filter.features_info )
                if filter.features_info(i).low_innovation_inlier == 1
                    featurectr=featurectr+1;
                    z(featurectr*2-1+(0:1)) = filter.features_info(i).z;
                    h(featurectr*2-1+(0:1)) = filter.features_info(i).h;
                    H(featurectr*2-1+(0:1), :)= filter.features_info(i).H;
                end
            end
            R = sparse(eye(length(z))*filter.sigmaCAM^2);
            filter.correctstates(h, z, H,R);
        end
        
        function rescue_hi_inliers( filter, cam )
            chi2inv_2_95 = 5.9915;
            % chi_099_2 = 9.2103;
            % predict points and compute H matrix
            predictPoints(filter, cam); % disable this when no low inlier
%             update is applied
            for i=1:length(filter.features_info)
                % do not use ~isempty(features_info(i).z), because for a
                % point of observation, z near the image border, 
                % predictPoints() may give empty prediction h and lostNo==1
                if ((filter.features_info(i).lostNo==0)&&(filter.features_info(i).low_innovation_inlier==0))
                    hi = filter.features_info(i).h; 
                    Si = filter.features_info(i).H*filter.p_k_k*filter.features_info(i).H';
                    nui = filter.features_info(i).z - hi;
                    if nui'/Si*nui<chi2inv_2_95
                        filter.features_info(i).high_innovation_inlier=1;
                    else
                        filter.features_info(i).high_innovation_inlier=0;
                    end
                end
            end
        end
        function SetAllHiInliers( filter)
            for i=1:length(filter.features_info)
                % do not use ~isempty(features_info(i).z), because for a
                % point of observation, z near the image border,
                % predictPoints() may give empty prediction h and lostNo==1
                if (filter.features_info(i).lostNo==0)
                    filter.features_info(i).high_innovation_inlier=1;
                end
            end
        end      

        function ekf_update_hi_inliers( filter )
            % mount vectors and matrices for the update
            allocsz=sum([filter.features_info.high_innovation_inlier]);
            if(allocsz==0)
%                 fprintf('No points in the states are high inliers!\n');
                return;
            end
            z=zeros(allocsz*2,1);
            h=z;
            H=zeros(allocsz*2, size(filter.p_k_k,1));
            featurectr=0;
            for i=1:length( filter.features_info)
                if filter.features_info(i).high_innovation_inlier == 1
                    featurectr=featurectr+1;
                    z(featurectr*2-1+(0:1)) = filter.features_info(i).z;
                    h(featurectr*2-1+(0:1)) = filter.features_info(i).h;
                    H(featurectr*2-1+(0:1), :)= filter.features_info(i).H;
                end
            end
            R = sparse(eye(length(z))*filter.sigmaCAM^2);
            filter.correctstates(h, z, H,R);
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
        % to be implemented later
        % for feature points are just lost, update the filter states with
        % their inverse depth triangulated from the Tracker measurements
        function ekf_update_rho(filter)
            zRho=zeros(1,length(filter.features_info));
            for smack=1:length(filter.features_info)
                if(filter.features_info(smack).lostNo==1 && isempty( filter.features_info(smack).z))
                    % just lost
                    zRho(smack)=filter.features_info(smack).invDepth;
                end
            end
            % update KF with measurements of Rho
        end
        % as we reduce the number of points in a group based on lostNo==1,
        % if we omit this, the lostNo may become 2 before it is used to reduce the point number in a group
              
        function   UpdateGrpPtsNo(filter, currentFrameNo)
            groupIds=[filter.groupPose.grpId];
            for i=1:length(filter.features_info)
                if(~isempty(filter.features_info(i).z))&&(filter.features_info(i).high_innovation_inlier==0)&&...
                        (filter.features_info(i).low_innovation_inlier==0 && (currentFrameNo-filter.features_info(i).initFrmNo)>15)
                    filter.features_info(i).lostNo=1; % speed up new feature incorporation
                end
                % lostNo and trackedNo may be increased by the tracker and
                % predictPoints()
                if(filter.features_info(i).lostNo==1)
                    % the feature point is just lost
                    grape=find(groupIds==filter.features_info(i).grpId,1);
                    filter.groupPose(grape).ptsNo=filter.groupPose(grape).ptsNo-1;
                    if(filter.groupPose(grape).ptsNo<0)
                        fprintf('Error in Renew, a group has negative number of points!');
                    end
                    filter.features_info(i).low_innovation_inlier = 0;
                    filter.features_info(i).high_innovation_inlier = 0;
                    filter.features_info(i).h = [];
                    filter.features_info(i).z = [];
                    filter.features_info(i).H = [];
                    filter.features_info(i).S = [];
                end  
                % do not touch z, h, S, H for other points as they are used
                % in update high inno inliers
            end
            % if the pts in the reference group frame is 0 or it has been
            % deleted, change the relay frame
            groupIds=[filter.groupPose.grpId]; % this groupIds may differ from the first groupIds
            tsar= find(groupIds==filter.refGrpFrmId,1);
            if(filter.groupPose(tsar).ptsNo==0)
                if((sum([filter.groupPose.ptsNo])~=0))
                    % find the group frame of the minimum variance and having
                    % points support
                    minCov=10000;
                    grpOrder=0;
                    for inch=1:length(filter.groupPose)
                        if(filter.groupPose(inch).ptsNo~=0)
                            spanner=filter.groupFrameSIP+(inch-1)*6;
                            cassette=max(diag(filter.p_k_k(spanner+(0:5),spanner+(0:5))));
                            if(minCov>cassette)
                                minCov=cassette;
                                grpOrder=inch;
                            end
                        end
                    end
                    % set its variance as epsilon so to fix the group frame
                    filter.refGrpFrmId=filter.groupPose(grpOrder).grpId;
                    % the following two methods are detrimental to the
                    % results
%                     spanner=filter.groupFrameSIP+(grpOrder-1)*6;
%                     rescale=100;
%                     filter.p_k_k(spanner+(0:5), spanner+(0:5))=filter.p_k_k(spanner+(0:5), spanner+(0:5))*rescale;
%                     filter.p_k_k(spanner+(0:5),:)=filter.p_k_k(spanner+(0:5),:)/rescale;
%                     filter.p_k_k(:, spanner+(0:5))=filter.p_k_k(:, spanner+(0:5))/rescale;                   
                    
%                     spanner=filter.groupFrameSIP+(grpOrder-1)*6;
%                     smallValue=0; %eps;
%                                     filter.p_k_k(spanner+(0:5),:)=smallValue;
%                                     filter.p_k_k(:, spanner+(0:5))=smallValue;
%                     filter.p_k_k(spanner+(0:5), spanner+(0:5))=smallValue;
                    % the 0 point support previous reference frame will be
                    % deleted once there is a new group coming
                end
            end            
        end
        % we must have sum(lostNo==0) == ptsNo in a group
        function CheckConstraints(filter)
            groupIds=[filter.groupPose.grpId];
            numPts=zeros(size(groupIds));
            for i=1:length(filter.features_info)              
                if(filter.features_info(i).lostNo==0)
                    % the feature point is not lost
                    grape=find(groupIds==filter.features_info(i).grpId,1);
                    numPts(grape)=numPts(grape)+1;
                end
            end
            assert(sum([filter.groupPose.ptsNo]-numPts)==0);
        end
        
        % update features_info, groupPose,
        % if there is a new group coming, remove the lost feature states and covariance, and remove the
        % lost groups and their covariance. Note this delayed
        % deletion of lost features and group frames does not cause trouble
        % in the filter since they have no measurement and remain stagnant
        % over the propagation phase.
        function   Renew( filter, newGrpPts)
            groupIds=[filter.groupPose.grpId];
            for i=1:length(filter.features_info)
                % lostNo and trackedNo were taken care by the tracker and
                % predictPoints()
                if(filter.features_info(i).lostNo==1)
                    % the feature point is just lost
                    grape=find(groupIds==filter.features_info(i).grpId,1);
                    filter.groupPose(grape).ptsNo=filter.groupPose(grape).ptsNo-1;
                    if(filter.groupPose(grape).ptsNo<0)
                        fprintf('Error in Renew, a group has negative number of points!');
                    end
                end
                
                filter.features_info(i).low_innovation_inlier = 0;
                filter.features_info(i).high_innovation_inlier = 0;
                filter.features_info(i).h = [];
                filter.features_info(i).z = [];
                filter.features_info(i).H = [];
                filter.features_info(i).S = [];
            end
            refFrameDeleted=~(length(filter.groupPose));
            if(~isempty(newGrpPts))
                %(1) delete lost feature points in the states
                deletion_list = zeros(length(filter.features_info),1);
                delctr=0;
                for i = 1:length(filter.features_info)
                    if (filter.features_info(i).lostNo>0)
                        delctr=delctr+1;
                        deletion_list(delctr) = i;
                    end
                end
                % Now, remove them from rear of the covariance and features_info vector
                filter.featSIP=filter.groupFrameSIP+length(filter.groupPose)*6;
                % the starting index of features' inverse depth in the
                % covariance matrix
                parToDelete=1;
                if delctr~=0
                    for i=delctr:-1:1
                        indexFromWichDeleteP=filter.featSIP+deletion_list(i)-1;
                        P_km1_km1_new = [filter.p_k_k(:,1:indexFromWichDeleteP-1) filter.p_k_k(:,indexFromWichDeleteP+parToDelete:end)];
                        filter.p_k_k = [P_km1_km1_new(1:indexFromWichDeleteP-1,:); P_km1_km1_new(indexFromWichDeleteP+parToDelete:end,:)];
                        if deletion_list(i)==1
                            filter.features_info = filter.features_info(2:end); continue;
                        end
                        if deletion_list(i)==length(filter.features_info)
                            filter.features_info = filter.features_info(1:end-1); continue;
                        end
                        if (deletion_list(i)~=length(filter.features_info))&&(deletion_list(i)~=1)
                            filter.features_info = [filter.features_info(1:deletion_list(i)-1) filter.features_info(deletion_list(i)+1:end)];
                        end
                    end
                    filter.p_k_k = sparse(filter.p_k_k);
                end
                %(2) remove groups of 0 support point and their covariance,
                % note (1) and (2) should be interchangable
                
                deletion_list = zeros(length(filter.groupPose),1);
                delctr=0;
                for i = 1:length(filter.groupPose)
                    if (filter.groupPose(i).ptsNo==0)
                        if(filter.groupPose(i).grpId==filter.refGrpFrmId)
                            refFrameDeleted=true;
                        end
                        delctr=delctr+1;
                        deletion_list(delctr) = i;
                    end
                end
                % Now, remove them from rear of the covariance and groupPose vector
                parToDelete=6;
                if delctr~=0
                    for i=delctr:-1:1
                        indexFromWichDeleteP=filter.groupFrameSIP+(deletion_list(i)-1)*6;
                        P_km1_km1_new = [filter.p_k_k(:,1:indexFromWichDeleteP-1) filter.p_k_k(:,indexFromWichDeleteP+parToDelete:end)];
                        filter.p_k_k = [P_km1_km1_new(1:indexFromWichDeleteP-1,:); P_km1_km1_new(indexFromWichDeleteP+parToDelete:end,:)];
                        if deletion_list(i)==1
                            filter.groupPose = filter.groupPose(2:end); continue;
                        end
                        if deletion_list(i)==length(filter.groupPose)
                            filter.groupPose = filter.groupPose(1:end-1); continue;
                        end
                        if (deletion_list(i)~=length(filter.groupPose))&&(deletion_list(i)~=1)
                            filter.groupPose = [filter.groupPose(1:deletion_list(i)-1) filter.groupPose(deletion_list(i)+1:end)];
                        end
                    end
                    filter.p_k_k=sparse(filter.p_k_k);
                end
                % add new feature points to the states and add the new
                % group frame to the groupPose, and update covariance
                grpPose.grpId=newGrpPts(1).grpId;
                grpPose.grpFrmNo=newGrpPts(1).initFrmNo;
                grpPose.pose=filter.rvqs0([7:10,1:3]);
                grpPose.ptsNo=length(newGrpPts);
                filter.groupPose=[filter.groupPose grpPose];% all the attributes form a column
                filter.features_info=[filter.features_info newGrpPts]; % the observations of an attribute form a row
                covDim=filter.groupFrameSIP-1+length(filter.groupPose)*6+length(filter.features_info);
                endIndex=filter.groupFrameSIP-1+length(filter.groupPose)*6-6; % the dim of previous states excluding the feature depths
                endIndex2=size(filter.p_k_k,2);
                Gmat=zeros(covDim, endIndex2+grpPose.ptsNo);
                Gmat(1:endIndex, 1: endIndex)=eye(endIndex);
                Gmat(endIndex+(1:3), filter.imuOrientSIP+(0:2))=eye(3); % for qs0 2 s
                Gmat(endIndex+(4:6), 7:9)=eye(3); % for Ts0 2s
                Gmat(endIndex+7:covDim, endIndex+1:end)=eye(-endIndex-6+covDim);
                Pk_pad=[filter.p_k_k, zeros(endIndex2, grpPose.ptsNo); 
                    zeros(grpPose.ptsNo, endIndex2), diag(ones(length(newGrpPts),1)*filter.std_rho^2)];
                filter.p_k_k= Gmat*Pk_pad*Gmat';
            end
            % if the pts in the reference group frame is 0 or it has been
            % deleted, change the relay frame, we put this procedure after the incorporation of
            % a new group because there may be no group of any support points
            groupIds=[filter.groupPose.grpId]; % this groupIds may differ from the first groupIds
            tsar= find(groupIds==filter.refGrpFrmId,1);
            if(refFrameDeleted||filter.groupPose(tsar).ptsNo==0)
                if(~isempty((filter.groupPose))&&(sum([filter.groupPose.ptsNo])~=0))
                    % find the group frame of the minimum variance and having
                    % points support
                    minCov=10000;
                    grpOrder=0;
                    for inch=1:length(filter.groupPose)
                        if(filter.groupPose(inch).ptsNo~=0)
                            spanner=filter.groupFrameSIP+(inch-1)*6;
                            cassette=max(diag(filter.p_k_k(spanner+(0:5),spanner+(0:5))));
                            if(minCov>cassette)
                                minCov=cassette;
                                grpOrder=inch;
                            end
                        end
                    end
                    % set its variance as epsilon so to fix the group frame
                    filter.refGrpFrmId=filter.groupPose(grpOrder).grpId;
%                     spanner=filter.groupFrameSIP+(grpOrder-1)*6;
%                     rescale=100;
%                     filter.p_k_k(spanner+(0:5), spanner+(0:5))=filter.p_k_k(spanner+(0:5), spanner+(0:5))*rescale;
%                     filter.p_k_k(spanner+(0:5),:)=filter.p_k_k(spanner+(0:5),:)/rescale;
%                     filter.p_k_k(:, spanner+(0:5))=filter.p_k_k(:, spanner+(0:5))/rescale;
%                     
                    % the 0 point support previous reference frame will be
                    % deleted once there is a new group coming
                end
            end
        end
        % return the feature depth and their columns in the feature table
        function ftJet=GetPtsDepth(filter)
            ftJet=[[filter.features_info.colFeatTable];[filter.features_info.invDepth]];
        end
        function SetFeaturesInfo(filter, trackedPoints)
            filter.features_info=trackedPoints;
        end
        function SaveToFile(filter, inillh_ant, preimutime, ffilres)
            Cen=llh2dcm_v000(ecef2geo_v000(filter.rvqs2e(1:3,1),0),[0;1]);
            vr_let=quat2dcm_v000(filter.rvqs2e(7:10));
            vr_c=rotro2eu('xyz',Cen*vr_let)*180/pi; % Cs2n
            xyz_ant=filter.rvqs2e(1:3)+quatrot_v000(filter.rvqs2e(7:10),filter.Tant2imu,0);
            if(~isempty(inillh_ant))
                vr_a=posdiff_v001(xyz_ant, inillh_ant);
                fwrite(ffilres,[preimutime;vr_a;filter.rvqs2e(4:6);vr_c;...
                    full(sqrt(diag(filter.p_k_k(7:15,7:15))))],'double');
            else
                fwrite(ffilres,[preimutime;ecef2geo_v000(xyz_ant,0);...
                    filter.rvqs2e(4:6);vr_c;full(sqrt(diag(filter.p_k_k(7:15,7:15))))],'double');
            end  
        end
        function SaveCamPoseandRqs02e(filter, imgTime, frmId)
                donkey=rotqr2eu('xyz', filter.camPose(1:4))*180/pi; % Cs2c
                mule=rotqr2eu('xyz', filter.rqs02e(4:7))*180/pi; % Cs2c
                fprintf(filter.camPoseFileHandle,...
                    '%d\t%15.8f\t%10.6f\t%10.6f\t%10.6f\t%8.5f\t%8.5f\t%8.5f\t%10.6f\t%10.6f\t%10.6f\t%8.5f\t%8.5f\t%8.5f\n',...
                    frmId, imgTime, donkey, filter.camPose(5:7),mule, filter.rqs02e(1:3));           
        end
        % to be implemented later
        % detect feature reappearance by predict its image coordinates
        % use normalized cross correlation to find the exact matching
        % point, put these points into the state vector and covariance,
        % measurement update with these points
        function LoopClosure(filter)
        end
        function mag=GetVelocityMag(filter)
            mag=norm(filter.rvqs2e(4:6),2);
        end        
    end
end