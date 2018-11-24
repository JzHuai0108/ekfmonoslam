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
function draw_camera( location, color )
% 
% DRAW_CAMERA: draw a camera in the current figure
%
% draw_camera( location, color )
%
% input
%   location: 3D position and orientation (via quaternion)

q = location(4:7);

camera_size = 0.2;

vertices = [ 0    1 -1  0 ;
             0    0  0  0 ;
            2 -1 -1 2  ]*camera_size;
        
% Rotate the camera
r = q2r(q);
vertices = r * vertices;

% Translate the camera
vertices = vertices + [location(1:3) location(1:3) location(1:3) location(1:3)];

% Draw the vertices
plot3(vertices(1,:), vertices(2,:), vertices(3,:), color, 'LineWidth', 2);