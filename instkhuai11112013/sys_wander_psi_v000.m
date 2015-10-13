%Psi-angle system model for wander mechanization under small angle
%assumption.
%Position errors are modeled as (Delta-R in meters) rather than
%((Delta-Llh) or (Delta-Theta)
%This is the wander angle implementation of the the psi-angle model described in Benson-1975

%Returns F and Q of discrete time model. The IMU error model must be
%specified here. (Note that in the previous implementations I had seperated
%nav and imu models in different m-files. I no longer follow that approach)

%This is a non-singular implementation
%For an example see: example_wander.m

function [STM Qd]=sys_wander_psi(Cen, h, Vn, qbn, acc, gyro, dt)

%%%Navigation error model
%1-3:Position error in meters
%4-6:vel errors
%7-9:Attitude errors in psi form

%N[da:dw]':Affect of IMU errors on nav states

Anav=zeros(9); %pos, vel, attitude
Nnav=zeros(9,6);

[Fc, wen_n, wie_n, g]=geoparam_v001(1, Cen(:,3), h, Vn);
Cbn=quat2dcm_v000(qbn);

%position errors
Anav(1:3,1:3)=-skew(wen_n);
Anav(1:3,4:6)=eye(3);

%Velocity Errors
Anav(4:5,1:2)=g*[-Fc(2,2) Fc(2,1);Fc(1,2) -Fc(1,1)];    %effect of horizontal position error on gravity
Anav(6,3)=2*g*Fc(1,1);  %effect of height on gravity (note that, drz is the negative of dh. Therefore, this is positive)
Anav(4:6,4:6)=-skew(wen_n+2*wie_n);
Anav(4:6,7:9)=skew(Cbn*acc);

Nnav(4:6,1:3)=Cbn; %%effect of accelerometer errors

%Attitude errors
Anav(7:9,7:9)=-skew(wen_n+ wie_n);
Nnav(7:9,4:6)=-Cbn; %Effect of gyroscope errors

%%%%Imu error model parameters
[Aimu_d, Qimu_d, Cimu, Rimu]=imu_err_model(acc, gyro, dt);

%%Convert continuous time Anav into discrete time
%I used van loan's method here to perform conversion so that I can obtain a
%correct Qnav_d. If you use a 1st order taylor series approximation, do not
%forget to correct the cross-correlation effects caused by arw/vrw. (In,
%the final implementation, it would be better to manually compute each
%element rather than using this method)
mx_a=dt*[-Anav,Nnav*Rimu*Nnav';zeros(9),Anav'];
mx_b = expm(mx_a);
Anav_d = mx_b(10:18,10:18)';
Qnav_d = Anav_d*mx_b(1:9,10:18);

%%Combine everything
nst_imu=size(Aimu_d);
STM=zeros(9+nst_imu);
Qd=zeros(9+nst_imu);

STM(1:9,1:9)=Anav_d;
STM(1:9,10:end)=Nnav*Cimu*dt;
STM(10:end,10:end)=Aimu_d;

Qd(1:9,1:9)=Qnav_d;
Qd(10:end,10:end)=Qimu_d;
Qd(1:9,10:end)=Nnav*Cimu*Qimu_d*dt/2;   %not necessary, but let's keep this
Qd(10:end,1:9)=Qd(1:9,10:end)';
end

function [Ad, Qd, C, R]=imu_err_model(acc, gyro, dt)
%%imu output model
%%x'=Ax+w, y=[acc;gyro]'+Cx+v, where <w,w>=Q, <v,v>=R
%%x(k+1)=Ad.x(k)+w(k) <w(k),w(k)>=Qd
%Note that the returning system model(Ad,Qd) is in discrete time, whereas the 
%bservation model (C,R) is in continuous time

%Imu error states
%1-3:Acc bias
%4-6:gyro bias
%7-9:Acc scale (ppt=part per thousand not million). ppm def can sometimes
%causes numerical problems
%10-12:Gyro scale

%load continuous time error model parameters
IMU_ERRDEF=imu_err_defs_v000(1);

%Discrete time model of imu error states for LN200
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

end