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

function uvu = undistort_fm( uvd, camera )
%
% Undistort image coordinates

% nPoints = size( uvd, 2 );
% uvu = zeros( 2, nPoints );
% for k = 1:nPoints;
%     uvu( :, k ) = undistor_a_point( uvd( :, k ), camera );
% end

Cx = camera.Cx;
Cy = camera.Cy;
k1 = camera.k1;
k2 = camera.k2;
dx = camera.dx;
dy = camera.dy;

xd = ( uvd(1,:) - Cx )*dx;
yd = ( uvd(2,:) - Cy )*dy;

rd = sqrt( xd.*xd + yd.*yd );

D = 1 + k1*rd.^2 + k2*rd.^4;
xu = xd.*D;
yu = yd.*D;

uvu = [ xu/dx + Cx; yu/dy + Cy ];