% system process model in local s0 frame with dcm formulation, Cs0e is the
% transformation from the sensor local s0 frame to e-frame
% rqs02e, r s0 in e, q s0 2e, rvqs0, rs in s0, vs in s0, qs0 2s
% acc, m/s^2, acc by imu in s frame, gyro, rad/s, angular rate by imu
% w.r.t i frame coordinated in s frame, dt, time interval for covariance update
% modelNo, imu model for bias and scale factor errors
function [STM Qd]=sys_local_dcm_bias(rqs02e, rvqs0, acc, gyro, dt, imutype, modelNo)
% the covariance corresponds to states, 
% rs in s0, v s in s0, q s02s, ba, bg
% for modelNo=3 and imutype=5, MEMS
% acc bias drift, random walk or random constant
% gyro bias drift, random walk or random constant

%system disturbance coefs
navStates=9;
N=zeros(navStates,6);
N(7:9,4:6)=eye(3); %attitude
Cs02s=quat2dcm_v000(rvqs0(7:10));
N(4:6,1:3)=Cs02s'; %velocity

%system matrix
A=zeros(navStates);
A(1:3, 4:6)=eye(3); % rs in s0
%Velocity
WIE_E=7292115e-11; 
wie2s0=quatrot_v000(rqs02e(4:7),[0;0;WIE_E], 1);
% A(4:6,1:3)=-skew(wie2s0)^2; % $\delta\bar{g}^{s0}(\mathbf{x}_s^{s0})$ is 
% noisy, but buried by low grade IMU noise
% A(4:6,4:6)=-2*skew(wie2s0); % for low grade IMU, w_ie^e is buried
A(4:6,7:9)=-Cs02s'*skew(acc);

%Attitude
A(7:9,7:9)=skew(-gyro);
Anav=A;
Nnav=N;
% X(k+1) = ffun[X(k),U(k),V(k)]
% X(k+1) = Ak*X(k)+Gk*Vk

%%%%Imu error model parameters
[Aimu_d, Qimu_d, Cimu, Rimu]=imu_err_model_v001(acc, gyro, dt, imutype, modelNo);
Aimu_d=Aimu_d(1:6,1:6);
Qimu_d=Qimu_d(1:6,1:6);
Cimu=Cimu(1:6,1:6);
%%%%Combine and discretize nav and imu models
% this discretization can also be accomplished by Loan's matrix exponential
% method, see sys_metric_phipsi_v000.m

Anav_d=eye(navStates)+dt*Anav;  %Use 1st order taylor series to discretize Anav
Qnav=Nnav*Rimu*Nnav';
Qnav_d=dt/2*(Anav_d*Qnav+Qnav*Anav_d');      %Use trapezoidal rule to discretize Rimu

STM=zeros(navStates+size(Aimu_d,1));
STM(1:navStates,1:navStates)=Anav_d;
STM(1:navStates,1+navStates:end)=Nnav*Cimu*dt;
STM(1+navStates:end,1+navStates:end)=Aimu_d;

Qd=zeros(navStates+size(Aimu_d,1));
Qd(1:navStates,1:navStates)=Qnav_d;
Qd(1+navStates:end,1+navStates:end)=Qimu_d;
Qd(1:navStates,1+navStates:end)=Nnav*Cimu*Qimu_d*dt/2;
Qd(1+navStates:end,1:navStates)=Qd(1:navStates,1+navStates:end)';
end