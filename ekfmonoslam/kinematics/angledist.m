function a = angledist(q1, q2)
% q1 [w, x, y, z]
% q2 [w, x, y, z]
q = quatmultiply(q1, quatconj(q2));
res = quat2axang(q);
a = res(4);
end
