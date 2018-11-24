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

function [x_k_k, p_k_k] = initialize_x_and_p()

% Initial velocity values
v_0 = 0;
std_v_0 = 0.025;
w_0 = 1e-15;
std_w_0 = 0.025; 

% Initial state vector and covariance matrix
x_k_k = [0 0 0 1 0 0 0 v_0 v_0 v_0 w_0 w_0 w_0]';
p_k_k = zeros(13,13);
p_k_k(1,1)=eps;
p_k_k(2,2)=eps;
p_k_k(3,3)=eps;
p_k_k(4,4)=eps;
p_k_k(5,5)=eps;
p_k_k(6,6)=eps;
p_k_k(7,7)=eps;
p_k_k(8,8)=std_v_0^2;
p_k_k(9,9)=std_v_0^2;
p_k_k(10,10)=std_v_0^2;
p_k_k(11,11)=std_w_0^2;
p_k_k(12,12)=std_w_0^2;
p_k_k(13,13)=std_w_0^2;