%-----------------------------------------------------------------------
% 1-point RANSAC EKF SLAM from a monocular sequence
%-----------------------------------------------------------------------

% Copyright (C) 2010 Javier Civera and J. M. M. Montiel
% Universidad de Zaragoza, Zaragoza, Spain.

% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation. Read http://www.gnu.org/copyleft/gpl.html for details

% If you use this code for academic work, please reference:
%   Javier Civera, Oscar G. Grasa, Andrew J. Davison, J. M. M. Montiel,
%   1-Point RANSAC for EKF Filtering: Application to Real-Time Structure from Motion and Visual Odometry,
%   to appear in Journal of Field Robotics, October 2010.

%-----------------------------------------------------------------------
% Authors:  Javier Civera -- jcivera@unizar.es 
%           J. M. M. Montiel -- josemari@unizar.es

% Robotics, Perception and Real Time Group
% Aragón Institute of Engineering Research (I3A)
% Universidad de Zaragoza, 50018, Zaragoza, Spain
% Date   :  May 2010
%-----------------------------------------------------------------------
function uv_rot=rotate_with_dist_fc(cam, uv_ini, R_c2c1, t_c2c1, n, d)
%
% trasfer the image of points throug a camera rotation and translation
%   the camera has radial distortion
% Input 
%   cam    - camera calibration
%   uv_ini - initial camera postion point
%   R_c2c1 - camera rotation matrix
% Output
%  uv_rot - points on the rotated image

uv_ini_und=undistort_fm(uv_ini',cam)';
uv_rot_und=inv(cam.K*(R_c2c1-(t_c2c1*n'/d))*inv(cam.K))*[uv_ini_und';ones(1,size(uv_ini_und,1))];
uv_rot_und=(uv_rot_und(1:2,:)./[uv_rot_und(3,:);uv_rot_und(3,:)])';
uv_rot=distort_fm(uv_rot_und',cam)';