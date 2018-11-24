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

function features_info = calculate_derivatives( x_k_km1, cam, features_info )

x_v = x_k_km1(1:13);
x_features = x_k_km1(14:end);

for i=1:length(features_info)
    
    if ~isempty(features_info(i).h)
        if strcmp(features_info(i).type, 'cartesian')
            y = x_features(1:3);
            x_features = x_features(4:end);
            features_info(i).H = sparse(calculate_Hi_cartesian( x_v, y, cam, i, features_info ));
        else
            y = x_features(1:6);
            x_features = x_features(7:end);
            features_info(i).H = sparse(calculate_Hi_inverse_depth( x_v, y, cam, i, features_info ));
        end
        
    else
        if strcmp(features_info(i).type, 'cartesian')
            x_features = x_features(4:end);
        end
        if strcmp(features_info(i).type, 'inversedepth')
            x_features = x_features(7:end);
        end
    end
    
end