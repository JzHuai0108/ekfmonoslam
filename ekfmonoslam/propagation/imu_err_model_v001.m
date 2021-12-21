% sys_metric_phipsi_v000.m recommends to split this function into a
% separate script file for reusability

% 1 without turn on bias estimates,first order GM bias and random walk scale factor errors
% 2 with random constant turn on bias estimates, first order GM bias and random walk scale factor errors
% 3 assumes bias and scale factors are random walks
% 4 assumes random constant bias and random constant scale factor errors
function [Ad, Qd, C, R]=imu_err_model_v001(acc, gyro, dt, imutype, modelNo)
%%imu output model
%%x'=Ax+w, y=[acc;gyro]'+Cx+v, where <w,w>=Q, <v,v>=R
%%x(k+1)=Ad.x(k)+w(k) <w(k),w(k)>=Qd

%NOTE THAT THE RETURNING SYSTEM MODEL(AD,QD) IS IN DISCRETE TIME, WHEREAS THE
%OBSERVATION MODEL (C,R) IS IN CONTINUOUS TIME
switch modelNo
    case 1
        %Note that the following lines are model definition dependent. Even if you use another imu model (for instance a model
        %withtout any scale factor error), do not modify these lines. Instead, just add a new "elseif" statement.
        %I know this looks wierd. However, a completely modular structure turns out to be impractible.
        %Imu error states
        %1-3:Acc bias, first order GM
        %4-6:gyro bias, first order GM
        %7-9:Acc scale (ppt=part per thousand not million). ppm def can
        %sometimes cause numerical problems, random walk
        %10-12:Gyro scale, random walk
        
        %load continuous time parameters
        IMU_ERRDEF=imu_err_defs_v000(imutype);
        
        %Discrete time model of imu error states
        acc_bias_a=exp(-dt/IMU_ERRDEF.acc_bias_Tc);
        acc_bias_w=IMU_ERRDEF.acc_bias_Tc*(IMU_ERRDEF.acc_bias_Q)/2*(1-acc_bias_a^2);  %(m/s^2)^2  (*s^2)
        gyro_bias_a=exp(-dt/IMU_ERRDEF.gyro_bias_Tc);
        gyro_bias_w=IMU_ERRDEF.gyro_bias_Tc*(IMU_ERRDEF.gyro_bias_Q)/2*(1-gyro_bias_a^2); %(rad/sec)^2 (*s^2)
        acc_scale_a=1;
        acc_scale_w=IMU_ERRDEF.acc_scale_Q*dt;     %ppt^2  (*s^2)
        gyro_scale_a=1;
        gyro_scale_w=IMU_ERRDEF.acc_scale_Q*dt;    %(ppt)^2  (*s^2)
        
        Ad=diag([acc_bias_a, acc_bias_a, acc_bias_a, gyro_bias_a, gyro_bias_a, gyro_bias_a, acc_scale_a, acc_scale_a, acc_scale_a, gyro_scale_a, gyro_scale_a, gyro_scale_a]);
        Qd=diag([acc_bias_w, acc_bias_w, acc_bias_w, gyro_bias_w, gyro_bias_w, gyro_bias_w, acc_scale_w, acc_scale_w, acc_scale_w, gyro_scale_w, gyro_scale_w, gyro_scale_w]);
        
        %Continuous time output model for the imu [acc;gyro]
        C=[eye(6), diag([acc(1),acc(2),acc(3),gyro(1),gyro(2),gyro(3)]/1000)];
        R=diag([IMU_ERRDEF.acc_vrw;IMU_ERRDEF.acc_vrw;IMU_ERRDEF.acc_vrw;IMU_ERRDEF.gyro_arw;IMU_ERRDEF.gyro_arw;IMU_ERRDEF.gyro_arw]);
    case 2 %mems 3dm gx3-35 sensor
        %Imu error states
        %1-3:Acc bias drift, first order GM process
        %4-6:gyro bias drift, first order GM process
        %7-9:Acc scale (ppt=part per thousand not million). ppm def can
        %sometimes cause numerical problems, random walk
        %10-12:Gyro scale factor, random walk
        %13-15: Acc turn on bias, random constant
        %16-18:gyro turn on bias, random constant
        %load continuous time parameters
        IMU_ERRDEF=imu_err_defs_v000(imutype);
        
        %Discrete time model of imu error states
        acc_bias_a=exp(-dt/IMU_ERRDEF.acc_bias_Tc);
        acc_bias_w=IMU_ERRDEF.acc_bias_Tc*(IMU_ERRDEF.acc_bias_Q)/2*(1-acc_bias_a^2);  %(m/s^2)^2  (*s^2)
        gyro_bias_a=exp(-dt/IMU_ERRDEF.gyro_bias_Tc);
        gyro_bias_w=IMU_ERRDEF.gyro_bias_Tc*(IMU_ERRDEF.gyro_bias_Q)/2*(1-gyro_bias_a^2); %(rad/sec)^2 (*s^2)
        acc_scale_a=1;
        acc_scale_w=IMU_ERRDEF.acc_scale_Q*dt;     %ppt^2  (*s^2)
        gyro_scale_a=1;
        gyro_scale_w=IMU_ERRDEF.acc_scale_Q*dt;    %(ppt)^2  (*s^2)
        
        Ad=diag([acc_bias_a, acc_bias_a, acc_bias_a, gyro_bias_a, gyro_bias_a, gyro_bias_a,acc_scale_a, acc_scale_a, acc_scale_a, gyro_scale_a, gyro_scale_a, gyro_scale_a,  ones(1,6)]);
        Qd=diag([acc_bias_w, acc_bias_w, acc_bias_w, gyro_bias_w, gyro_bias_w, gyro_bias_w,acc_scale_w, acc_scale_w, acc_scale_w, gyro_scale_w, gyro_scale_w, gyro_scale_w, zeros(1,6)]);
        
        %Continuous time output model for the imu [acc;gyro]
        C=[eye(6), diag([acc(1),acc(2),acc(3),gyro(1),gyro(2),gyro(3)]/1000), eye(6)];
        R=diag([IMU_ERRDEF.acc_vrw;IMU_ERRDEF.acc_vrw;IMU_ERRDEF.acc_vrw;IMU_ERRDEF.gyro_arw;IMU_ERRDEF.gyro_arw;IMU_ERRDEF.gyro_arw]);
        
    case 3
        %  Imu error states
        %1-3:Acc bias drift, random walk
        %4-6:gyro bias drift, random walk
        %7-9:Acc scale (ppt=part per thousand not million). ppm def can
        %sometimes cause numerical problems, random walk
        %10-12:Gyro scale factor, random walk
        %load continuous time parameters
        IMU_ERRDEF=imu_err_defs_v000(imutype);
        
        %Discrete time model of imu error states
        acc_bias_a=1;
        acc_bias_w=IMU_ERRDEF.acc_bias_Q*dt;  %(m/s^2)^2  (*s^2)
        gyro_bias_a=1;
        gyro_bias_w=IMU_ERRDEF.gyro_bias_Q*dt; %(rad/sec)^2 (*s^2)
        acc_scale_a=1;
        acc_scale_w=IMU_ERRDEF.acc_scale_Q*dt;     %ppt^2  (*s^2)
        gyro_scale_a=1;
        gyro_scale_w=IMU_ERRDEF.acc_scale_Q*dt;    %(ppt)^2  (*s^2)
        
        Ad=diag([acc_bias_a, acc_bias_a, acc_bias_a, gyro_bias_a, gyro_bias_a, gyro_bias_a, acc_scale_a, acc_scale_a, acc_scale_a, gyro_scale_a, gyro_scale_a, gyro_scale_a]);
        Qd=diag([acc_bias_w, acc_bias_w, acc_bias_w, gyro_bias_w, gyro_bias_w, gyro_bias_w, acc_scale_w, acc_scale_w, acc_scale_w, gyro_scale_w, gyro_scale_w, gyro_scale_w]);
        
        %Continuous time output model for the imu [acc;gyro]
        C=[eye(6), diag([acc(1),acc(2),acc(3),gyro(1),gyro(2),gyro(3)]/1000)];
        R=diag([IMU_ERRDEF.acc_vrw;IMU_ERRDEF.acc_vrw;IMU_ERRDEF.acc_vrw;IMU_ERRDEF.gyro_arw;IMU_ERRDEF.gyro_arw;IMU_ERRDEF.gyro_arw]);
        
        
    case 4 % constant bias and constant scale factor for high end IMUs
        %Imu error states
        %1-3:Acc bias drift, random constant
        %4-6:gyro bias drift, random constant
        %7-9:Acc scale (ppt=part per thousand not million). ppm def can
        %sometimes cause numerical problems, random constant
        %10-12:Gyro scale factor, random constant
        
        IMU_ERRDEF=imu_err_defs_v000(imutype);        
        %Discrete time model of imu error states     
        Ad=sparse(eye(12));        
        Qd=zeros(12);         
        %Continuous time output model for the imu [acc;gyro]
        C=[eye(6), diag([acc(1),acc(2),acc(3),gyro(1),gyro(2),gyro(3)]/1000)];
        R=diag([IMU_ERRDEF.acc_vrw;IMU_ERRDEF.acc_vrw;IMU_ERRDEF.acc_vrw;IMU_ERRDEF.gyro_arw;IMU_ERRDEF.gyro_arw;IMU_ERRDEF.gyro_arw]);
    
    case 5
        %  Imu error states
        %1-3:Acc bias drift, random walk
        %4-6:gyro bias drift, random constant
        %7-9:Acc scale (ppt=part per thousand not million). ppm def can
        %sometimes cause numerical problems, random constant
        %10-12:Gyro scale factor, random constant
        %load continuous time parameters
        IMU_ERRDEF=imu_err_defs_v000(imutype);
        
        %Discrete time model of imu error states
        acc_bias_a=1;
        acc_bias_w=IMU_ERRDEF.acc_bias_Q*dt;  %(m/s^2)^2  (*s^2)
        gyro_bias_a=1;
        gyro_bias_w=0; %(rad/sec)^2 (*s^2)
        acc_scale_a=1;
        acc_scale_w=0;     %ppt^2  (*s^2)
        gyro_scale_a=1;
        gyro_scale_w=0;    %(ppt)^2  (*s^2)
        
        Ad=diag([acc_bias_a, acc_bias_a, acc_bias_a, gyro_bias_a, gyro_bias_a, gyro_bias_a, acc_scale_a, acc_scale_a, acc_scale_a, gyro_scale_a, gyro_scale_a, gyro_scale_a]);
        Qd=diag([acc_bias_w, acc_bias_w, acc_bias_w, gyro_bias_w, gyro_bias_w, gyro_bias_w, acc_scale_w, acc_scale_w, acc_scale_w, gyro_scale_w, gyro_scale_w, gyro_scale_w]);
        
        %Continuous time output model for the imu [acc;gyro]
        C=[eye(6), diag([acc(1),acc(2),acc(3),gyro(1),gyro(2),gyro(3)]/1000)];
        R=diag([IMU_ERRDEF.acc_vrw;IMU_ERRDEF.acc_vrw;IMU_ERRDEF.acc_vrw;IMU_ERRDEF.gyro_arw;IMU_ERRDEF.gyro_arw;IMU_ERRDEF.gyro_arw]);
        
        
    otherwise
        disp('This is an undefined error model.');
        
end
