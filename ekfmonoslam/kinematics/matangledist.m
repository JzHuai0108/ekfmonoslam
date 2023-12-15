function a = matangledist(m1, m2)
% m1 4x4 transformation matrix or 3x3 rotation matrix
res = m1(1:3, 1:3) * transpose(m2(1:3, 1:3));
aa = rotm2axang(res);
a = aa(4);
end
