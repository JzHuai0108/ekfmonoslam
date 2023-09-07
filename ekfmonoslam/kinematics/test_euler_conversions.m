function test_euler_conversions()
% Thankfully, ekfmonoslam, instk, oxford robotcar dataset devkit, MATLAB
% are using the same Euler conventions.
% As an aside, ROS tf transform, Eigen, scipy.spatial.transform.Rotation
% also use the same Euler convention as can be seen from 
% https://github.com/JzHuai0108/vio_common/blob/master/python/test/test_euler.py

eul = rand(3, 1) * 2 - 1; % rpy
format longg;
R = euler2dcm_v000(eul);

R2 = eul2rotm_oxford(eul(1), eul(2), eul(3));
assert(norm(R - R2) < 3e-7);

eul2 = dcm2euler_v000(R);
assert(norm(eul - eul2) < 3e-7);
eul3 = rotro2eu('xyz', R);
assert(norm(eul2 - eul3) < 3e-7);

eul4 = Cbn2att(R);
assert(norm(eul3 - eul4) < 3e-7);

quat = quaternion(R, 'rotmat', 'point');
eulZYX = quat2eul(quat, 'zyx'); % a row vector
assert(norm(fliplr(eulZYX) - eul3') < 3e-7);

Rp = quat2rotm(quat);
assert(norm(Rp - R) < 3e-7);
end
