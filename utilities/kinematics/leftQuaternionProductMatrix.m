function pl = leftQuaternionProductMatrix(p)
% p 4x1 w,x,y,z
% p * q = leftQuaternionProductMatrix(p) * q
% eq 18 in Sola Quaternion kinematics for the error-state Kalman filter
w = p(1);
x = p(2);
y = p(3);
z = p(4);
pl = [w, -x, -y, -z; 
    x, w, -z, y;
    y, z, w, -x;
    z, -y, x, w];
end
