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

function features_info = add_feature_to_info_vector( uv, im_k, X_RES, features_info, step, newFeature )

half_patch_size_when_initialized = 20;
half_patch_size_when_matching = 6;
new_position = length(features_info)+1;

features_info(new_position).patch_when_initialized =...
    im_k(uv(2,1)-half_patch_size_when_initialized:uv(2,1)+half_patch_size_when_initialized,...
    uv(1,1)-half_patch_size_when_initialized:uv(1,1)+half_patch_size_when_initialized);
features_info(new_position).patch_when_matching = zeros( 2*half_patch_size_when_matching+1, 2*half_patch_size_when_matching+1 );
features_info(new_position).r_wc_when_initialized=X_RES(1:3);
features_info(new_position).R_wc_when_initialized=q2r(X_RES(4:7));
features_info(new_position).uv_when_initialized=uv';
features_info(new_position).half_patch_size_when_initialized=half_patch_size_when_initialized;
features_info(new_position).half_patch_size_when_matching=half_patch_size_when_matching;
features_info(new_position).times_predicted=0;
features_info(new_position).times_measured=0;
features_info(new_position).init_frame = step;
features_info(new_position).init_measurement = uv;
features_info(new_position).type = 'inversedepth';
features_info(new_position).yi = newFeature;
features_info(new_position).individually_compatible = 0;
features_info(new_position).low_innovation_inlier = 0;
features_info(new_position).high_innovation_inlier = 0;
features_info(new_position).z = [];
features_info(new_position).h = [];
features_info(new_position).H = [];
features_info(new_position).S = [];
features_info(new_position).state_size = 6;
features_info(new_position).measurement_size = 2;
features_info(new_position).R = eye(features_info(new_position).measurement_size);