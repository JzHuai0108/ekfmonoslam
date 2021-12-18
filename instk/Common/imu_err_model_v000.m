%%This script converts application specific IMU models into a standard
%%format that can be used by sys_model scripts.

%%Do not change the existing models!!!!!.
%%Whenever you need an additional model, just define a new imutype and add
%%it to this script.
%%SEE also the notes in "imu_err_defs_v000"

function [Ad, Qd, C, R]=imu_err_model(acc, gyro, dt, imutype)
%%imu output model
%%x'=Ax+w, y=[acc;gyro]'+Cx+v, where <w,w>=Q, <v,v>=R
%%x(k+1)=Ad.x(k)+w(k) <w(k),w(k)>=Qd
%Note that the returning system model(Ad,Qd) is in discrete time, whereas the 
%bservation model (C,R) is in continuous time

%load continuous time parameters
IMU_ERRDEF=imu_err_defs_v000(imutype);

if imutype==0 %dummy to obtain nav system model from sys_XXX scripts
    Ad=[];
    Qd=[];
    C=[];
    R=zeros(6);
elseif ismember(imutype, [1,2,3])
    %Imu error states
    %1-3:Acc bias
    %4-6:gyro bias
    %7-9:Acc scale (ppt=part per thousand not million). ppm def can sometimes
    %causes numerical problems
    %10-12:Gyro scale
    
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
elseif ismember(imutype, 4)
    %Imu error states
    %1-3:Acc bias
    %4-6:gyro bias
    
    acc_bias_a=exp(dt*IMU_ERRDEF.Acc.A);
    gyro_bias_a=exp(dt*IMU_ERRDEF.Gyro.A);
    acc_bias_w=IMU_ERRDEF.Acc.Q/(-IMU_ERRDEF.Acc.A*2)*(1-acc_bias_a^2);  %(m/s^2)^2  (*s^2)
    gyro_bias_w=IMU_ERRDEF.Gyro.Q/(-IMU_ERRDEF.Gyro.A*2)*(1-gyro_bias_a^2);  %(m/s^2)^2  (*s^2)
    
    Ad=diag([ones(1,3)*acc_bias_a, ones(1,3)*gyro_bias_a]);
    Qd=diag([ones(1,3)*acc_bias_w, ones(1,3)*gyro_bias_w]);
    
    C=eye(6);
    R=diag([ones(1,3)*IMU_ERRDEF.Acc.R, ones(1,3)*IMU_ERRDEF.Gyro.R]);
else
    disp('sys_metric_phipsi: IMUTYPE is not defined');
end
end