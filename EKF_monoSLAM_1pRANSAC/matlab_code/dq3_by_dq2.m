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

function dq3_by_dq2RES=dq3_by_dq2(q1_in)

q1.R=q1_in(1);
q1.X=q1_in(2);
q1.Y=q1_in(3);
q1.Z=q1_in(4);

dq3_by_dq2RES = ...
   [q1.R, -q1.X, -q1.Y, -q1.Z,
    q1.X,  q1.R,  q1.Z, -q1.Y,
    q1.Y, -q1.Z,  q1.R,  q1.X,
    q1.Z,  q1.Y, -q1.X,  q1.R];

return
