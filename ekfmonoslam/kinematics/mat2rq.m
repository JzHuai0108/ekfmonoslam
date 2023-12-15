function rq = mat2rq(m)
% input: a 4x4 transformation matrix
% output: rq [x, y, z, qx, qy, qz, qw]
q = rotm2quat(m(1:3, 1:3));
rq =[m(1:3, 4)', q(4), q(1:3)];
end
