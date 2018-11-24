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

function features_info = set_as_most_supported_hypothesis( features_info, positions_li_inliers_id, positions_li_inliers_euc )

j_id = 1;
j_euc = 1;

for i=1:length(features_info)
    
    if ~isempty(features_info(i).z)
        if strcmp(features_info(i).type, 'cartesian')
            if positions_li_inliers_euc(j_euc)
                features_info(i).low_innovation_inlier = 1;
            else
                features_info(i).low_innovation_inlier = 0;
            end
            j_euc = j_euc + 1;
        end
        if strcmp(features_info(i).type, 'inversedepth')
            if positions_li_inliers_id(j_id)
                features_info(i).low_innovation_inlier = 1;
            else
                features_info(i).low_innovation_inlier = 0;
            end
            j_id = j_id + 1;
        end
    end
    
end