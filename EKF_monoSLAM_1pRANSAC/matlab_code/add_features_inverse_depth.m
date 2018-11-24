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

function [ X_RES, P_RES, newFeature ] = add_features_inverse_depth( uvd, X, P, cam, std_pxl, initial_rho, std_rho )

nNewFeat = size( uvd, 2 );

if nNewFeat == 0
    X_RES = X;
    P_RES = P;
    return;

else

    % Camera state
    Xv = X(1:13);

    X_RES = X;
    P_RES = P;

    for i = 1:nNewFeat
        newFeature = hinv( uvd(:,i), Xv, cam, initial_rho );
        X_RES = [ X_RES; newFeature ];
        P_RES = add_a_feature_covariance_inverse_depth( P_RES, uvd(:,i), Xv, std_pxl, std_rho, cam );
    end

end