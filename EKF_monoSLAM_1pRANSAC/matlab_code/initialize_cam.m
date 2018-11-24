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

function cam = initialize_cam()

d =     0.0112;
nRows = 240;
nCols = 320;
Cx =    1.7945 / d;
Cy =    1.4433 / d;
k1=     6.333e-2;
k2=     1.390e-2;
f =     2.1735;

cam.k1 =    k1;
cam.k2 =    k2;
cam.nRows = nRows;
cam.nCols = nCols;
cam.Cx =    Cx;
cam.Cy =    Cy;
cam.f =     f;
cam.dx =    d;
cam.dy =    d;
cam.model = 'two_distortion_parameters';

cam.K =     sparse( [ f/d   0     Cx;
                0  f/d    Cy;
                0    0     1] );