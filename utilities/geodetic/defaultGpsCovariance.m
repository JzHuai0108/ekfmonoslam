function R = defaultGpsCovariance(quality, hdop)
% The default variances for solutions of different quality refers to 
% https://github.com/ros-drivers/nmea_navsat_driver/blob/master/src/libnmea_navsat_driver/driver.py
if nargin < 2
    hdop = 2;
end
% epe corresponds to quality number 0,1,2,3,4,5,6,7,8,9
default_epe = [1000000, 4.0, 0.1, 1000000, 0.02, 4.0, 1000000, 1000000, 1000000, 3.0];
if (quality >= 0 && quality <= 9)
    e = default_epe(quality + 1);
    R = diag([e * hdop; e * hdop; 4 * e * hdop].^2);
else
    e = default_epe(1);
    R = diag([e * hdop; e * hdop; 4 * e * hdop].^2);
end
end
