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

function [ cartesian ] = inversedepth2cartesian( inverse_depth )

rw = inverse_depth(1:3,:);
theta = inverse_depth(4,:);
phi = inverse_depth(5,:);
rho = inverse_depth(6,:);
% m = mfunc(theta, phi,1);
cphi = cos(phi);
m = [cphi.*sin(theta); -sin(phi); cphi.*cos(theta)];
cartesian(1,:) = rw(1) + (1./rho).*m(1,:);
cartesian(2,:) = rw(2) + (1./rho).*m(2,:);
cartesian(3,:) = rw(3) + (1./rho).*m(3,:);