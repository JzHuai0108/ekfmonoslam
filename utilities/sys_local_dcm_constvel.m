% system process model in local s0 frame with dcm formulation, Cs0e is the
% transformation from the sensor local s0 frame to e-frame
% written by Huai
% rqs02e, r s0 in e, q s0 2e, rvqs0, rs in s0, vs in s0, qs0 2s
% acc, m/s^2, acc by imu in s frame, gyro, rad/s, angular rate by imu
% w.r.t i frame coordinated in s frame, dt, time interval for covariance update
% modelNo, imu model for bias and scale factor errors
function [STM Qd]=sys_local_dcm_constvel(rqs02e, rvqs0, acc, gyro, dt, imutype, velNoiseStd)
% the covariance corresponds to states, 
% rs in s0, v s in s0, q s02s, ba, bg
% this function assumes velocity noise, disregarding accelerometer input

navStates=9;
%system matrix
A=zeros(navStates);
A(1:3, 4:6)=eye(3); % rs in s0
%Velocity
WIE_E=7292115e-11; 
wie2s0=quatrot_v000(rqs02e(4:7),[0;0;WIE_E], 1);
% A(4:6,1:3)=-skew(wie2s0)^2; % $\delta\bar{g}^{s0}(\mathbf{x}_s^{s0})$ is 
% noisy, but buried by low grade IMU noise
% A(4:6,4:6)=-2*skew(wie2s0); % for low grade IMU, w_ie^e is buried
% A(4:6,7:9)=-Cs02s'*skew(acc);

%Attitude
A(7:9,7:9)=skew(-gyro);

%%Combine everything
nst_imu=6;
STM=zeros(9+nst_imu);
Qd=zeros(9+nst_imu);
Gk=zeros(9+nst_imu,6); % velocity noise and gyro ARW
STM(1:9,1:9)=A;
STM(7:9,13:15)=eye(3);
STM=eye(9+nst_imu)+STM*dt; % FIRST ORDER TAYLOR APPROXIMATION

Gk( 4: 6, 1: 3)	= eye(3);		% noise driving velocity 
Gk( 7: 9, 4: 6)	= eye(3);		% ARW
IMU_ERRDEF=imu_err_defs_v000(imutype);
R=diag([velNoiseStd^2;velNoiseStd^2; (velNoiseStd*2)^2;IMU_ERRDEF.gyro_arw;IMU_ERRDEF.gyro_arw;IMU_ERRDEF.gyro_arw]);
Qd=Gk*R*dt*Gk';
end