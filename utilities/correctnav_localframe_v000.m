%Corrects the navigation states using the Kalman generated error values

%dx=corrections, in order of rqs02e, rvqs0
%delta_x(k)=K(delta_y(k)-h(delta_x(k|k-1)))=K*delta_y(k)
% the complete state is defined as delta(ecef x,y,z of sensor at epoch t0,
% q s0 2 e frame(s0 frame of sensor is tied to epoch t0 of the imu)),
% the delta(pos s in s0 frame, vx, vy, vz in s0 frame of sensor, q s0 to s frame)
% accelerometer and gyro bias, and
% acc scale and gyro scale factor
function [rqs02e, rvqs0]=correctnav_localframe_v000(rqs02e, rvqs0, dx)

%Correct the position
rqs02e(1:3)=rqs02e(1:3)-dx(1:3);
qet=rvec2quat_v000(dx(4:6));
rqs02e(4:7)=quatmult_v000(qet, rqs02e(4:7));

rvqs0(1:3)=rvqs0(1:3)-dx(7:9);
%velocity correction
rvqs0(4:6)=rvqs0(4:6)-dx(10:12);
qst=rvec2quat_v000(dx(13:15));
rvqs0(7:10)=quatmult_v000(qst, rvqs0(7:10));
%     %Method II:
%     %Method I involves some trigonometric functions. Here is a
%     %algebraic correction method (See: Chung 1996-Eq:16) (Special thanks to Kaygisiz)
%     R=[-q(2) -q(3) -q(4);q(1) q(4) -q(3);q(-4) q(1) q(2);q(3) q(-2) q(1)];
%     dq=-0.5*R*dx(7:9);
%     att_new=att_new-dq;
%     att_new=att_new/norm(att_new);  %optional
end
