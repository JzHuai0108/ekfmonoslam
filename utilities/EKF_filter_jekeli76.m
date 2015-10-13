
classdef EKF_filter_jekeli76 < handle
    properties (Hidden)
        type = 'ekf';
        tag  = 'EKF_crude';  % ID tag
        covDim=3; % the dimension of covariance matrix, Rx RY RZ in earth centered
        % n-frame, vn, RPY, acc bias, gyro bias, acc scale and gyro scale 
    end
    % The following properties can be set only by class methods
    properties (SetAccess = private)   
        dt; % sampling interval of IMU, unit sec
        % the translation from antenna to mems frame, i.e., the antenna's
        % position in the mems frame.
        xi;
        vi;
        bias=0;  
        p_k_k; % the covariance of the entire state vector  
        psdq; % process noise PSD
    end
    methods
        function filter = EKF_filter_jekeli76(initX, initV, deltat, procnoiseq)   
            filter.xi=initX;
            filter.vi=initV;
            filter.dt=deltat;
            filter.psdq=procnoiseq;
            filter.p_k_k=eye(filter.covDim);      
        end
        %===============================================================================================
        %-- State transition function
        % propogate state with accelerometer and gyro input at time k-1 to state at k,
        % i.e., X(k|k-1), from state at k-1. U1 contains IMU measurement  
        % U1 is the measured acceleration
        function ffun_state(filter, U1)
            
            ahat=U1-filter.bias;
            vi1=filter.vi+ahat*filter.dt;            
           filter.xi=filter.xi+(vi1+filter.vi)/2*filter.dt;
           filter.vi=vi1;            
        end
        function ffun_covariance(filter )
            %propagate the covariance corresponds to states, rs in e, vs in e, q s2e,
            deltat=filter.dt;
            STM=[1, 0, deltat; deltat, 1, .5*deltat^2; 0,0, 1];
            Qd=filter.psdq*[deltat, .5*deltat^2, 0; .5*deltat^2, deltat^3/3, 0; 0, 0,0];
            filter.p_k_k=STM*filter.p_k_k*STM'+Qd;  % the covariance of the navigation states and imu error terms
        end
        %==============================================================================================
        %update both the covariance and state
        function correctstates(filter, predict,measure, H,R)
            p_km1_k=filter.p_k_k;
            inno= predict-measure;
            %Kalman
            K=p_km1_k*H'/(H*p_km1_k*H'+R);
            deltaX=K*inno;
            % update covariance
            filter.p_k_k=(eye(filter.covDim)-K*H)*p_km1_k*(eye(filter.covDim)-K*H)'+K*R*K';
            % compute updated state
            filter.xi=filter.xi-deltaX(1);
            filter.vi=filter.vi-deltaX(2);
            filter.bias = filter.bias + deltaX(3);
        end
        function SaveToFile(filter, preimutime, ffilres)
            %Write result to the files
            fprintf(ffilres,'%6.2f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\n', ...
                [preimutime;filter.vi;filter.xi;filter.bias;sqrt(diag(filter.p_k_k))]);         
        end
    end
end