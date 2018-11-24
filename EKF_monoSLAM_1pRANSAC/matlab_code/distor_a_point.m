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

function uvd = distor_a_point( uvu, camera )

  Cx = camera.Cx;
  Cy = camera.Cy;
  k1 = camera.k1;
  k2 = camera.k2;
  dx = camera.dx;
  dy = camera.dy;
  
  
  xu=(uvu(1)-Cx)*dx;
  yu=(uvu(2)-Cy)*dy;
  
  ru=sqrt(xu*xu+yu*yu);
  rd=ru/(1+k1*ru^2+k2*ru^4);
  
  % 10 iterations are enough to distort the usual monoSLAM camera image.
  % Other cameras, run check_distortion.m
  n_iterations = 20;
  
  for k=1:n_iterations
      f=rd+k1*rd^3+k2*rd^5-ru;
      f_p=1+3*k1*rd^2+5*k2*rd^4;
      rd=rd -f/f_p;
  end
    
  D=1+k1*rd^2+k2*rd^4;
  xd=xu/D;
  yd=yu/D;
  
  uvd=[xd/dx+Cx; yd/dy+Cy];
