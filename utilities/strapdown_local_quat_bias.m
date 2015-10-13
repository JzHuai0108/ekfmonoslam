% input: rvqs0: r s in s0, v s in s0, q s0 2 s,
% rqs02e: rs0 in e, q s0 2 e, a, acceleration by IMU, w, angular rate by
% gyro, dt, sample interval from the last measurement to the current
% measurement.
% gwomegaw, $g^w$, gravity in local world frame, 3-vector, and $w_{iw}^w$, 3-vector,
% if all of them set 0, this amount to integration in inertial frame

function rvqs0_new=strapdown_local_quat_bias(rvqs0, rqs02e, a, w, dt, gwomegaw)
rvqs0_new=zeros(size(rvqs0));
gs0=gwomegaw (1:3); % gravity in s0 frame
wie2s0=gwomegaw (4:6);
%Update attitude
%% method (1) second order integration
qe=rvec2quat_v000(wie2s0*dt);
vr_a=quatmult_v000(rvqs0(7:10),qe);
qb=rvec2quat_v000(-w*dt);
rvqs0_new(7:10)=quatmult_v000(qb,vr_a);
%% method (2) Runge-Kutta 4th order integration, empirically, this sometimes
% gives worse result than second order integration
% wie2s=quatrot_v000(rvqs0(7:10),wie2s0,0);
% omega=zeros(4,2);
% omega(1,2)=dt;
% omega(2:4,1)=w-wie2s;
% omega(2:4,2)=lastw-wie2s;
% qs2s0=rotationRK4( omega, [rvqs0(7); -rvqs0(8:10)]);
% rvqs0_new(7:10)=[qs2s0(1); -qs2s0(2:4)];

%% better velocity and position integration than first order rectanglar rule
%Update Velocity
vel_inc1=(quatrot_v000(rvqs0(7:10),a*dt,1)+quatrot_v000(rvqs0_new(7:10),a*dt,1))/2;
vel_inc2=(gs0+2*cross(rvqs0(4:6),wie2s0))*dt;
rvqs0_new(4:6)=rvqs0(4:6)+vel_inc1+vel_inc2;
%Update_pos
rvqs0_new(1:3)=rvqs0(1:3)+(rvqs0(4:6)+rvqs0_new(4:6))*dt/2;
