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

function [ X_km1_k, P_km1_k ] = predict_state_and_covariance( X_k, P_k, type, SD_A_component_filter, SD_alpha_component_filter )


delta_t = 1;

% camera motion prediction
Xv_km1_k = fv( X_k(1:13,:), delta_t, type, SD_A_component_filter, SD_alpha_component_filter  );

% features prediction
X_km1_k = [ Xv_km1_k; X_k( 14:end,: ) ];

% state transition equation derivatives
F = sparse( dfv_by_dxv( X_k(1:13,:),zeros(6,1),delta_t, type ) );% civera
% model
% F=eye(12); % huai's model
% state noise
linear_acceleration_noise_covariance = (SD_A_component_filter*delta_t)^2;
angular_acceleration_noise_covariance = (SD_alpha_component_filter*delta_t)^2;
Pn = sparse (diag( [linear_acceleration_noise_covariance linear_acceleration_noise_covariance linear_acceleration_noise_covariance...
        angular_acceleration_noise_covariance angular_acceleration_noise_covariance angular_acceleration_noise_covariance] ) );

Q = func_Q( X_k(1:13,:), zeros(6,1), Pn, delta_t, type); % civera model
% G= [eye(3)*delta_t, zeros(3);zeros(3), eye(3)*delta_t; eye(3), zeros(3);zeros(3), eye(3)];
% Q=G*Pn*G'; % huai's model

size_P_k = size(P_k,1);
% huai's parameters
% P_km1_k = [ F*P_k(1:12,1:12)*F' + Q         F*P_k(1:12,13:size_P_k);
%             P_k(13:size_P_k,1:12)*F'        P_k(13:size_P_k,13:size_P_k)];
% civera's model
P_km1_k = [ F*P_k(1:13,1:13)*F' + Q         F*P_k(1:13,14:size_P_k);
            P_k(14:size_P_k,1:13)*F'        P_k(14:size_P_k,14:size_P_k)];