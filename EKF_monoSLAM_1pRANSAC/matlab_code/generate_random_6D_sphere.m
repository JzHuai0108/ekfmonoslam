% generate random 6d [x,y,z, theta, phi,rho] vectors in each column, in
% total 1000 columns and 6 rows
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

nPointsRand = 1000;
X = rand( 1, nPointsRand )-0.5;
Y = rand( 1, nPointsRand )-0.5;
Z = rand( 1, nPointsRand )-0.5;
theta = rand( 1, nPointsRand )-0.5;
phi = rand( 1, nPointsRand )-0.5;
lambda = rand( 1, nPointsRand )-0.5;
for i = 1:nPointsRand
    a = [X(i) Y(i) Z(i) theta(i) phi(i) lambda(i)];
    a = a/norm(a)*sqrt(12.59158724374398);
    X(i) = a(1); Y(i) = a(2); Z(i) = a(3);
    theta(i) = a(4); phi(i) = a(5); lambda(i) = a(6);
end
randSphere6D = [ X; Y; Z; theta; phi; lambda ];