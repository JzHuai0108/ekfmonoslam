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

function d=dq_by_deuler(euler_angles);

phi = euler_angles(1);
theta = euler_angles(2);
psi = euler_angles(3);

d = [(0.5)*(-sin(phi/2)+cos(phi/2))     (0.5)*(-sin(theta/2)+cos(theta/2))      (0.5)*(-sin(psi/2)+cos(psi/2)); 
     (0.5)*(+cos(phi/2)+sin(phi/2))     (0.5)*(-sin(theta/2)-cos(theta/2))      (0.5)*(-sin(psi/2)-cos(psi/2)); 
     (0.5)*(-sin(phi/2)+cos(phi/2))     (0.5)*(+cos(theta/2)-sin(theta/2))      (0.5)*(-sin(psi/2)+cos(psi/2)); 
     (0.5)*(-sin(phi/2)-cos(phi/2))     (0.5)*(-sin(theta/2)-cos(theta/2))      (0.5)*(+cos(psi/2)+sin(psi/2))];