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

function [X_km1_km1_new,P_km1_km1_new] = delete_a_feature( X_km1_km1,P_km1_km1,featToDelete,features_info )


if strcmp(features_info(featToDelete).type,'cartesian')
    parToDelete = 3;
else
    parToDelete = 6;
end

indexFromWichDelete = 14;

for i=1:featToDelete-1
    if strcmp(features_info(i).type, 'inversedepth')
        indexFromWichDelete = indexFromWichDelete + 6;
    end
    if strcmp(features_info(i).type, 'cartesian')
        indexFromWichDelete = indexFromWichDelete + 3;
    end
end

X_km1_km1_new = [X_km1_km1(1:indexFromWichDelete-1); X_km1_km1(indexFromWichDelete+parToDelete:end)];
    

P_km1_km1_new = [P_km1_km1(:,1:indexFromWichDelete-1) P_km1_km1(:,indexFromWichDelete+parToDelete:end)];
P_km1_km1_new = [P_km1_km1_new(1:indexFromWichDelete-1,:); P_km1_km1_new(indexFromWichDelete+parToDelete:end,:)];