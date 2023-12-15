function pose = rq2mat(d)
    % input: d = [x, y, z, qx, qy, qz, qw]
    % output: transformation matrix 4x4
    pose = eye(4);
    q = quaternion(d(7), d(4), d(5), d(6));
    pose(1:3, 1:3) = quat2rotm(q);
    pose(1:3, 4) = [d(1), d(2), d(3)]';
end