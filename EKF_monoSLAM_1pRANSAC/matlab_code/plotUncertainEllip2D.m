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

function plotUncertainEllip2D(C,nu,chi2,color,linewidth)

hold on;

th=0:2*pi/100:2*pi;

x=[cos(th') sin(th')]'*sqrt(chi2);

if(min(eig(C))<0)
    C=eye(2);
    color=[0 0 0];
    fprintf('NPSD matrix, a black false ellipse has been plot\n');
end
K=chol(C)';

nPoints=size(th,2);
y=K*x+[ones(1,nPoints)*nu(1);ones(1,nPoints)*nu(2)];

h = plot(y(1,:),y(2,:));

set(h,'Color',color,'LineWidth',linewidth);