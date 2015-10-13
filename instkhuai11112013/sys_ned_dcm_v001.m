%System model for a geodesic frame navigation based on PHI formulation. 
%att must be either Cbn or qbn.
%Position errors are modelled as delta_Llh (not in metric)

%This script is the same as sys_ned_dcm_v000 except the fact that imu error
%models are incorporated within this script. (_v000 is coded to be called
%by mdl_ned_dcm. I no longer follow that approach)

function [STM Qd]=sys_ned_dcm(Llh, vel_n, att, acc, gyro, dt, imutype)
% position   (1-3)
% velocity   (4-6)
% attitude   (7-9)

if (size(att,1)==1 || size(att,2)==1) %att is qbn
    Cbn=quat2dcm_v000(att);
else %att is dcm
    Cbn=att;
end

%%%%Navigation Error parameters
%system disturbance coefs
Nnav=zeros(9,6);
Nnav(7:9,4:6)=-Cbn; %attitude
Nnav(4:6,1:3)=Cbn; %velocity

%system matrix
Anav=zeros(9);
[Rn, Re, ~, sL, cL, WIE_E]=geoparam_v000(Llh);
tL=sL/cL;
Rn_h=Rn+Llh(3);
Re_h=Re+Llh(3);
acc_n=Cbn*acc;

Anav(1,3)=-vel_n(1)/(Rn_h)^2; 
Anav(1,4)=1/(Rn_h);

Anav(2,1)=vel_n(2)*tL/(Re_h)/cL;
Anav(2,3)=-vel_n(2)/(Re_h)/(Re_h)/cL;
Anav(2,5)=1/(Re_h)/cL;

Anav(3,6)=-1;

Anav(4,1)=-(2*WIE_E*cL*vel_n(2)+vel_n(2)*vel_n(2)/(Re_h)/cL/cL);
Anav(4,3)=(-vel_n(1)*vel_n(3)/(Rn_h)/(Rn_h)+vel_n(2)*vel_n(2)*tL/(Re_h)/(Re_h));
Anav(4,4)=vel_n(3)/(Rn_h);
Anav(4,5)=-(2*WIE_E*sL+2*vel_n(2)*tL/(Re_h));
Anav(4,6)=vel_n(1)/(Rn_h);
Anav(4,8)=-acc_n(3);
Anav(4,9)=acc_n(2);

Anav(5,1)=(2*WIE_E*cL*vel_n(1)-2*WIE_E*sL*vel_n(3)+vel_n(2)*vel_n(1)/(Re_h)/cL/cL);
Anav(5,3)=-(vel_n(2)*vel_n(3)+vel_n(2)*vel_n(1)*tL)/(Re_h)/(Re_h);
Anav(5,4)=(2*WIE_E*sL+vel_n(2)*tL/(Re_h));
Anav(5,5)=(vel_n(3)+vel_n(1)*tL)/(Re_h);
Anav(5,6)=(2*WIE_E*cL+vel_n(2)/(Re_h));
Anav(5,7)=acc_n(3);
Anav(5,9)=-acc_n(1);

Anav(6,1)=2*WIE_E*sL*vel_n(2);
Anav(6,3)=(vel_n(1)*vel_n(1)/(Rn_h)/(Rn_h)+vel_n(2)*vel_n(2)/(Re_h)/(Re_h));
Anav(6,4)=-2*vel_n(1)/(Rn_h);
Anav(6,5)=-2*(WIE_E*cL+vel_n(2)/(Re_h));
Anav(6,7)=-acc_n(2);
Anav(6,8)=acc_n(1);


Anav(7,1)=-WIE_E*sL;
Anav(7,3)=-vel_n(2)/(Re_h)/(Re_h);
Anav(7,5)=1/(Re_h);
Anav(7,8)=-(WIE_E*sL+vel_n(2)*tL/(Re_h));
Anav(7,9)=vel_n(1)/(Rn_h);

Anav(8,3)=vel_n(1)/(Rn_h)/(Rn_h);
Anav(8,4)=-1/(Rn_h);
Anav(8,7)=(WIE_E*sL+vel_n(2)*tL/(Re_h));
Anav(8,9)=(WIE_E*cL+vel_n(2)/(Re_h));

Anav(9,1)=-(WIE_E*cL+vel_n(2)/(Re_h)/cL/cL);
Anav(9,3)=vel_n(2)*tL/(Re_h)/(Re_h);
Anav(9,5)=-tL/(Re_h);
Anav(9,7)=-vel_n(1)/(Rn_h);
Anav(9,8)=-(WIE_E*cL+vel_n(2)/(Re_h));


%%%%Imu error model parameters
[Aimu_d, Qimu_d, Cimu, Rimu]=imu_err_model(acc, gyro, dt, imutype);

%%%%Combine and discretize nav and imu models
Anav_d=eye(9)+dt*Anav;  %Use 1st order taylor series to discretize Anav
Qnav=Nnav*Rimu*Nnav';
Qnav_d=dt/2*(Anav_d*Qnav+Qnav*Anav_d');      %Use trapezoidal rule to discretize Rimu

STM=zeros(21,21);
STM(1:9,1:9)=Anav_d;
STM(1:9,10:end)=Nnav*Cimu*dt;
STM(10:end,10:end)=Aimu_d;

Qd=zeros(21);
Qd(1:9,1:9)=Qnav_d;
Qd(10:end,10:end)=Qimu_d;
Qd(1:9,10:end)=Nnav*Cimu*Qimu_d*dt/2;
Qd(10:end,1:9)=Qd(1:9,10:end)';

end


function [Ad, Qd, C, R]=imu_err_model(acc, gyro, dt, imutype)
%%imu output model
%%x'=Ax+w, y=[acc;gyro]'+Cx+v, where <w,w>=Q, <v,v>=R
%%x(k+1)=Ad.x(k)+w(k) <w(k),w(k)>=Qd

%NOTE THAT THE RETURNING SYSTEM MODEL(AD,QD) IS IN DISCRETE TIME, WHEREAS THE 
%OBSERVATION MODEL (C,R) IS IN CONTINUOUS TIME

%Imu error states
%1-3:Acc bias
%4-6:gyro bias
%7-9:Acc scale (ppt=part per thousand not million). ppm def can sometimes
%causes numerical problems
%10-12:Gyro scale

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

end