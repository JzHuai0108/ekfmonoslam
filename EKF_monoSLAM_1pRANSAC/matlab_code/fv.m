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

function X_k_km1=fv(X_k_k,delta_t, type, std_a, std_alpha)

rW =X_k_k(1:3,1);
qWR=X_k_k(4:7,1);
vW =X_k_k(8:10,1);
wW =X_k_k(11:13,1);

if strcmp(type,'constant_orientation')
    wW = [0 0 0]';
    X_k_km1=[rW+vW*delta_t;
        qWR;
        vW;
        wW];
end

if strcmp(type,'constant_position')
    vW = [0 0 0]';
    X_k_km1=[rW;
        reshape(qprod(qWR,v2q(wW*delta_t)),4,1);
        vW;
        wW];
end

if strcmp(type,'constant_position_and_orientation')
    vW = [0 0 0]';
    wW = [0 0 0]';
    X_k_km1=[rW;
        qWR;
        vW;
        wW];
end

if strcmp(type,'constant_position_and_orientation_location_noise')
    vW = [0 0 0]';
    wW = [0 0 0]';
    X_k_km1=[rW;
        qWR;
        vW;
        wW];
end

if strcmp(type,'constant_velocity')
    X_k_km1=[rW+vW*delta_t;
%         quatmult_v000(qWR, rvec2quat_v000(wW*delta_t));
        reshape(qprod(qWR,v2q(wW*delta_t)),4,1);
        vW;
        wW];
%     % verfication
%     tempq1=quatmult_v000(qWR, rvec2quat_v000(wW*delta_t));
%     tempq2=reshape(qprod(qWR,v2q(wW*delta_t)),4,1);
%     temp=tempq2-tempq1;
%     disp(temp');
end