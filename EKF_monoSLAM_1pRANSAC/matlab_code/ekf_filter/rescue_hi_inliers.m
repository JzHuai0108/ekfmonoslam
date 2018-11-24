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

function features_info = rescue_hi_inliers( filter, features_info, cam )

chi2inv_2_95 = 5.9915;
% chi_099_2 = 9.2103;

features_info = predict_camera_measurements( filter.x_k_k, cam, features_info );
features_info = calculate_derivatives( filter.x_k_k, cam, features_info );

for i=1:length(features_info)
    
    if (features_info(i).individually_compatible==1)&&(features_info(i).low_innovation_inlier==0)
        hi = features_info(i).h';
        Si = features_info(i).H*filter.p_k_k*features_info(i).H';
        nui = features_info(i).z - hi;
        if nui'*inv(Si)*nui<chi2inv_2_95
            features_info(i).high_innovation_inlier=1;
        else
            features_info(i).high_innovation_inlier=0;
        end
    end
end