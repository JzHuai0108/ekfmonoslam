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

function uvu = undistor_a_point( uvd, camera )

Cx = camera.Cx;
Cy = camera.Cy;
k1 = camera.k1;
k2 = camera.k2;
dx = camera.dx;
dy = camera.dy;

ud = uvd(1);
vd = uvd(2);
rd = sqrt( ( dx*(ud-Cx) )^2 + (dy*(vd-Cy) )^2 );

uvu = [ Cx + ( ud - Cx )*( 1 + k1*rd^2 + k2*rd^4 ); Cy + ( vd - Cy )*( 1 + k1*rd^2 + k2*rd^4 ) ];