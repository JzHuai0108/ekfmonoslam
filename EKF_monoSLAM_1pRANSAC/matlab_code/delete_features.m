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

function [filter, features_info ] = delete_features( filter, features_info )

% First, list all features that should be deleted
deletion_list = [];
for i = 1:length(features_info)
    if ( features_info(i).times_measured < 0.5*features_info(i).times_predicted ) &&...
            ( features_info(i).times_predicted > 5 )
        deletion_list = [deletion_list; i];
    end
end

% Now, remove them from the state, covariance and features_info vector
if size(deletion_list,1)~=0

    for i=size(deletion_list,1):-1:1
            [ x_k_k, p_k_k ] = delete_a_feature( get_x_k_k(filter),...
                get_p_k_k(filter), deletion_list(i), features_info );
            filter = set_x_k_k(filter,x_k_k);
            filter = set_p_k_k(filter,p_k_k);

        if deletion_list(i)==1
            features_info = features_info(2:end);
        end

        if deletion_list(i)==size(features_info,2)
            features_info = features_info(1:end-1);
        end

        if (deletion_list(i)~=size(features_info,2))&&(deletion_list(i)~=1)
            features_info = [features_info(1:deletion_list(i)-1) features_info(deletion_list(i)+1:end)];
        end
    end
end