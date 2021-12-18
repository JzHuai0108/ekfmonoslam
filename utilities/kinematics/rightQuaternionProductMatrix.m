function qr = rightQuaternionProductMatrix(q)
% q 4x1 w,x,y,z
% p * q = rightQuaternionProductMatrix(q) * p
% eq 18 in Sola Quaternion kinematics for the error-state Kalman filter
w = q(1);
x = q(2);
y = q(3);
z = q(4);
qr = [w, -x, -y, -z; 
    x, w, z, -y;
    y, -z, w, x;
    z, y, -x, w];
end
