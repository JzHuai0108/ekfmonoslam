% system process model in local s0 frame with dcm formulation, Cs0e is the
% transformation between the sensor local s0 frame and e-frame
% written by Huai
% rqs02e, r s0 in e, q s0 2e, rvqs0, rs in s0, vs in s0, qs0 2s
% acc, m/s^2, acc by imu in s frame, gyro, rad/s, angular rate by imu
% w.r.t i frame coordinated in s frame, dt, time interval for covariance update
function [STM Qd]=sys_local_dcm_v000(rqs02e, rvqs0, acc, gyro, dt, imutype, modelNo)
% the covariance corresponds to states, rs0 in e, q s02e,
% rs in s0, v s in s0, q s02s, ba, bg, sa, sg, qs2c, Ts in c,
% for each group frame, qs02si and Tsi in s0, for each
% point, its inverse depth, here we are only interested up to sa and sg
% for modelNo=3 and imutype=5, MEMS
% acc bias drift, random walk
% gyro bias drift, random walk
% acc scale factor, random walk
% gyro scale factor, random walk

% for modelNo=1 and imutype=5 the acc and gyro turn on bias are removed
R=6317000; % earth average radius
%system disturbance coefs
N=zeros(15,6);
N(13:15,4:6)=eye(3); %attitude
Cs02s=quat2dcm_v000(rvqs0(7:10));
N(10:12,1:3)=Cs02s'; %velocity

%system matrix
A=zeros(15);
% 0s for r s0 in e and q s0 to e
A(7:9, 10:12)=eye(3); % rs in s0
%Velocity
xyz_imu=rqs02e(1:3)+quatrot_v000(rqs02e(4:7),rvqs0(1:3),0);
unit=xyz_imu/norm(xyz_imu,2);
Llh=ecef2geo_v000(xyz_imu,0);
[Rn, Re, g, sL, cL, WIE_E]=geoparam_v000(Llh);
geCoeff=-g*R^2/norm(xyz_imu,2)^3*(eye(3)-3*(unit*unit'))-skew([0;0;WIE_E])*skew([0;0;WIE_E]);
Cs02e=quat2dcm_v000(rqs02e(4:7));
A(10:12,1:3)=Cs02e'*geCoeff;
A(10:12,7:9)=Cs02e'*geCoeff*Cs02e;

Cne=pos2Cne_v000(Llh(1), Llh(2));
ge=Cne(:,3)*g-[xyz_imu(1:2);0]*WIE_E^2; % gravitation in e frame
A(10:12,4:6)=Cs02e'*(-skew(ge)+geCoeff*skew(quatrot_v000(rqs02e(4:7),rvqs0(1:3),0))+...
    skew(skew([0;0;WIE_E])*skew([0;0;WIE_E])*xyz_imu)-2*skew(quatrot_v000(rqs02e(4:7),rvqs0(4:6),0))*skew([0;0;WIE_E]));

A(10:12, 10:12)=-2*skew(quatrot_v000(rqs02e(4:7),[0;0;WIE_E], 1));
A(10:12, 13:15)=-Cs02s'*skew(acc);
%Attitude
A(13:15,4:6)=Cs02s*Cs02e'*skew([0;0;WIE_E]);
A(13:15,13:15)=skew(-gyro);
Anav=A;
Nnav=N;
% X(k+1) = ffun[X(k),U(k),V(k)]
% X(k+1) = Ak*X(k)+Gk*Vk

%%%%Imu error model parameters
[Aimu_d, Qimu_d, Cimu, Rimu]=imu_err_model_v001(acc, gyro, dt, imutype, modelNo);

%%%%Combine and discretize nav and imu models
% this discretization can also be accomplished by Loan's matrix exponential
% method, see sys_metric_phipsi_v000.m
navStates=15;
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