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

function uv_u = hu( yi, cam)

u0 = cam.Cx;
v0 = cam.Cy;
f  = cam.f;% f in mm, ku in number of pixels per mm
ku = 1/cam.dx;
kv = 1/cam.dy;

uv_u = zeros( 2, size( yi, 2 ) );

for i = 1:size( yi, 2 )
    uv_u( 1, i ) = u0 + (yi(1,i)/yi(3,i))*f*ku; 
    uv_u( 2, i ) = v0 + (yi(2,i)/yi(3,i))*f*kv;
end