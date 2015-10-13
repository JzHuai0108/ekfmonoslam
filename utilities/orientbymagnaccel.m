function euler=orientbymagnaccel(accel, mag, V)
% given accelerometer reading(m/s^2) in sensor frame and magnetometer 
% reading( Gauss) in sensor frame, determine the euler angles corresponding
% to Rs2n(sensor frame to n frame) such that Rn2s=R1(roll)R2(pitch)R3(yaw),
% return [roll; pitch; yaw]. yaw is relative to magnetic north
% optionally, hard iron offset estimate V may be provided for accuracy
% This is a rough implementation of the method in Implementing a Tilt-
% Compensated eCompass using Accelerometer and Magnetometer Sensors
% http://cache.freescale.com/files/sensors/doc/app_note/AN4248.pdf 
accel=-accel;% follow exactly "Implementing a Tilt-...."
phi=atan2(accel(2), accel(3));
theta=atan2(-accel(1), accel(2)*sin(phi)+accel(3)*cos(phi));
if(nargin==1)
    euler=[phi; theta; 0];
    return;
end
if(nargin==3)
    Bpmv=mag-V;
else
    Bpmv=mag;
end
psi=atan2(Bpmv(3)*sin(phi)-Bpmv(2)*cos(phi), Bpmv(1)*cos(theta)+...
    Bpmv(2)*sin(theta)*sin(phi)+Bpmv(3)*sin(theta)*cos(phi));
euler=[phi; theta; psi];
end