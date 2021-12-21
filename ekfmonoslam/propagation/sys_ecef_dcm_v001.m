% system process model in ecef frame with dcm formulation, Cse is the 
% transformation between the sensor and e-frame, with the assumption that
% the sensor is also the body frame, we have Cse=Cbe.
%implement this function as sys_llh_phipsi_v000
function [STM Qd]=sys_ecef_dcm_v001(xyz_imu, Cse, acc, gyro, dt, imutype, modelNo)
% for modelNo=5 and imutype=5
% position   (1-3)
% velocity   (4-6)
% attitude   (7-9) 
% acc bias drift
% gyro bias drift
% acc scale factor
% gyro scale factor
% acc turn on bias
% gyro turn on bias
% for modelNo=1 and imutype=5 the acc and gyro turn on bias are removed
WIE_E=7292115e-11;  %Earth rotation rate
g = 9.7803267714; % nominal g on earth surface
R=6317000; % earth average radius
%system disturbance coefs
N=zeros(9,6);
N(7:9,4:6)=-Cse; %attitude
N(4:6,1:3)=Cse; %velocity

%system matrix
A=zeros(9);
acc_e=Cse*acc;

%Position
A(1,4)=1;
A(2,5)=1;
A(3,6)=1;

%Velocity
unit=xyz_imu/norm(xyz_imu,2);
A(4:6,1:3)=-g*R^2/norm(xyz_imu,2)^3*(eye(3)-3*(unit*unit'))-skew([0;0;WIE_E])*skew([0;0;WIE_E]);
A(4:6,4:6)=-2*skew([0;0;WIE_E]);
A(4:6,7:9)=skew(acc_e);

%Attitude
A(7:9,7:9)=-skew([0;0;WIE_E]);
Anav=A;
Nnav=N;
% X(k+1) = ffun[X(k),U(k),V(k)]
% X(k+1) = Ak*X(k)+Gk*Vk

%%%%Imu error model parameters
[Aimu_d, Qimu_d, Cimu, Rimu]=imu_err_model_v001(acc, gyro, dt, imutype, modelNo);

%%%%Combine and discretize nav and imu models
% this discretization can also be accomplished by Loan's matrix exponential
% method, see sys_metric_phipsi_v000.m
Anav_d=eye(9)+dt*Anav;  %Use 1st order taylor series to discretize Anav
Qnav=Nnav*Rimu*Nnav';
Qnav_d=dt/2*(Anav_d*Qnav+Qnav*Anav_d');      %Use trapezoidal rule to discretize Rimu

STM=zeros(9+size(Aimu_d,1));
STM(1:9,1:9)=Anav_d;
STM(1:9,10:end)=Nnav*Cimu*dt;
STM(10:end,10:end)=Aimu_d;

Qd=zeros(9+size(Aimu_d,1));
Qd(1:9,1:9)=Qnav_d;
Qd(10:end,10:end)=Qimu_d;
Qd(1:9,10:end)=Nnav*Cimu*Qimu_d*dt/2;
Qd(10:end,1:9)=Qd(1:9,10:end)';
end

