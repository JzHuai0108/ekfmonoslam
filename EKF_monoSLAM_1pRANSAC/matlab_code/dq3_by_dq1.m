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

function dq3_by_dq1RES=dq3_by_dq1(q2_in)

 q2.R=q2_in(1); 
 q2.X=q2_in(2);
 q2.Y=q2_in(3);
 q2.Z=q2_in(4);
 
 dq3_by_dq1RES=[q2.R, -q2.X, -q2.Y, -q2.Z,
                q2.X,  q2.R, -q2.Z,  q2.Y,
		        q2.Y,  q2.Z,  q2.R, -q2.X,
		        q2.Z, -q2.Y,  q2.X,  q2.R];
 return
