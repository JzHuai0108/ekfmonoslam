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

function uv_c1=rotate_with_dist_fc_c1c2(cam, uv_c2, R_c1c2, t_c1c2, n, d)
%
% trasfer the image of points throug a camera rotation and translation
%   the camera has radial distortion
% Input 
%   cam    - camera calibration
%   uv_ini - initial camera postion point
%   R_c1c2 - camera rotation matrix
% Output
%  uv_rot - points on the rotated image

K = cam.K;

uv_c2_und=undistort_fm(uv_c2',cam)';
uv_c1_und=K*(R_c1c2-(t_c1c2*n'/d))*inv(K)*[uv_c2_und';ones(1,size(uv_c2_und,1))];
uv_c1_und=(uv_c1_und(1:2,:)./[uv_c1_und(3,:);uv_c1_und(3,:)])';
uv_c1=distort_fm(uv_c1_und',cam)';