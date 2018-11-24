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

function [ filter, features_info ] = inversedepth_2_cartesian( filter, features_info )

linearity_index_threshold = 0.1;
X = get_x_k_k(filter);
P = get_p_k_k(filter);

for i=1:length(features_info)

    if strcmp(features_info(i).type, 'inversedepth')
        initialPositionOfFeature = 14;
        for j=1:i-1
            if strcmp(features_info(j).type, 'cartesian')initialPositionOfFeature = initialPositionOfFeature + 3; end
            if strcmp(features_info(j).type, 'inversedepth') initialPositionOfFeature = initialPositionOfFeature + 6; end
        end
        std_rho = sqrt(P(initialPositionOfFeature + 5,initialPositionOfFeature + 5));
        rho = X(initialPositionOfFeature + 5);
        std_d = std_rho/(rho^2);
        % camera 2 to point distance
        theta = X(initialPositionOfFeature + 3);
        phi = X(initialPositionOfFeature + 4);
        mi = m(theta,phi);
        x_c1 = X(initialPositionOfFeature:initialPositionOfFeature+2);
        x_c2 = X(1:3);
        p = inversedepth2cartesian( X(initialPositionOfFeature:initialPositionOfFeature+5) );
        d_c2p = norm(p - x_c2);
        % alpha (parallax)
        cos_alpha = ((p-x_c1)'*(p-x_c2))/(norm(p-x_c1)*norm(p-x_c2));

        % Linearity index
        linearity_index = 4*std_d*cos_alpha/d_c2p;

        if linearity_index<linearity_index_threshold
            size_X_old = size(X,1);
            % State Vector
            X = [X(1:initialPositionOfFeature-1); p; X(initialPositionOfFeature+6:end)];
            % Covariance
            dm_dtheta = [cos(phi)*cos(theta)  0   -cos(phi)*sin(theta)]';
            dm_dphi = [-sin(phi)*sin(theta)  -cos(phi)   -sin(phi)*cos(theta)]';
            J = [ eye(3) (1/rho)*dm_dtheta (1/rho)*dm_dphi -mi/(rho^2) ];
            J_all = sparse( [eye(initialPositionOfFeature-1) zeros(initialPositionOfFeature-1,6) zeros(initialPositionOfFeature - 1, size_X_old - initialPositionOfFeature - 5);...
                zeros(3, initialPositionOfFeature-1) J zeros(3, size_X_old - initialPositionOfFeature - 5);...
                zeros(size_X_old - initialPositionOfFeature - 5, initialPositionOfFeature - 1) zeros(size_X_old - initialPositionOfFeature - 5, 6) eye(size_X_old - initialPositionOfFeature - 5)] );
            P = J_all*P*J_all';
            filter = set_x_k_k(filter,X);
            filter = set_p_k_k(filter,P);
            features_info(i).type = 'cartesian';
            return; % for the moment, only convert one feature per step (sorry!)
        end
    end
end