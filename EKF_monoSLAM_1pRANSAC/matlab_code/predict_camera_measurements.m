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

function features_info = predict_camera_measurements( x_k_k, cam, features_info )

% Pinhole model
t_wc = x_k_k(1:3);
r_wc = q2r(x_k_k(4:7));
features = x_k_k( 14:end );

for i = 1:length(features_info)
    
    if strcmp(features_info(i).type, 'cartesian')
        yi = features( 1:3 );
        features = features( 4:end );
        hi = hi_cartesian( yi, t_wc, r_wc, cam, features_info(i) );
        if (~isempty(hi))
            features_info(i).h = hi';
        end
    end
    
    if strcmp(features_info(i).type, 'inversedepth')
        yi = features( 1:6 );
        features = features( 7:end );
        hi = hi_inverse_depth( yi, t_wc, r_wc, cam, features_info(i) );
        if (~isempty(hi))
            features_info(i).h = hi';
        end
    end
    
end