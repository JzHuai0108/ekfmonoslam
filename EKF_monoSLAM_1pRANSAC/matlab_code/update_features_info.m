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

function features_info = update_features_info( features_info )

% reset features_info parameters
for i=1:length(features_info)
    if ~isempty(features_info(i).h)
        features_info(i).times_predicted = features_info(i).times_predicted + 1;
    end
    if (features_info(i).low_innovation_inlier||features_info(i).high_innovation_inlier)
        features_info(i).times_measured = features_info(i).times_measured + 1;
    end
    if(isempty(features_info(i).h)&&isempty(features_info(i).z))
        features_info(i).dormancy=features_info(i).dormancy+1;
    elseif(~isempty(features_info(i).z))
        features_info(i).dormancy=0;
    end
    features_info(i).individually_compatible = 0;
    features_info(i).low_innovation_inlier = 0;
    features_info(i).high_innovation_inlier = 0;
    features_info(i).h = [];
    features_info(i).z = [];
    features_info(i).H = [];
    features_info(i).S = [];
end