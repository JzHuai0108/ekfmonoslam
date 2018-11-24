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

function [state_vector_pattern, z_id, z_euc] = generate_state_vector_pattern( features_info, x )

state_vector_pattern = zeros(length(x),4);
position = 14;
z_id = [];
z_euc = [];

for i=1:length( features_info )
    
    if strcmp(features_info(i).type, 'inversedepth')
        if ~isempty(features_info(i).z)
            state_vector_pattern(position:position+2,1) = ones(3,1);
            state_vector_pattern(position+3:position+4,2) = ones(2,1);
            state_vector_pattern(position+5,3) = 1;
            z_id = [z_id [features_info(i).z(1); features_info(i).z(2)]];
        end
        position = position + 6;
    end
    if strcmp(features_info(i).type, 'cartesian')
        if ~isempty(features_info(i).z)
            state_vector_pattern(position:position+2,4) = ones(3,1);
            z_euc = [z_euc [features_info(i).z(1); features_info(i).z(2)]];
        end
        position = position + 3;
    end
    
end