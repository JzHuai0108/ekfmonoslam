function d = matrdist(m1, m2)
% m1 is a 4x4 transformation matrix.
d = norm(m1(1:3, 4) - m2(1:3, 4), 2); 
end
