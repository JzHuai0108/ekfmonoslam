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

function J_undistor=jacob_undistor_fm(camera,uvd)
%
% Jacobian of the undistortion of the image coordinates
%  presented in
%  Real-Time 3D SLAM with Wide-Angle Vision, 
%      Andrew J. Davison, Yolanda Gonzalez Cid and Nobuyuki Kita, IAV 2004.
% input
%    camera   -  camera calibration parameters
%    uvd       - distorted image points in pixels
% output
%    J_dunistor -  Jacobian

Cx=camera.Cx;
Cy=camera.Cy;
k1=camera.k1;
k2=camera.k2;
dx=camera.dx;
dy=camera.dy;
  
ud=uvd(1);
vd=uvd(2);
xd=(uvd(1)-Cx)*dx;
yd=(uvd(2)-Cy)*dy;
  
rd2=xd*xd+yd*yd;
rd4=rd2*rd2;
     
uu_ud=(1+k1*rd2+k2*rd4)+(ud-Cx)*(k1+2*k2*rd2)*(2*(ud-Cx)*dx*dx);
vu_vd=(1+k1*rd2+k2*rd4)+(vd-Cy)*(k1+2*k2*rd2)*(2*(vd-Cy)*dy*dy);
    
uu_vd=(ud-Cx)*(k1+2*k2*rd2)*(2*(vd-Cy)*dy*dy);
vu_ud=(vd-Cy)*(k1+2*k2*rd2)*(2*(ud-Cx)*dx*dx);
     
J_undistor=[uu_ud uu_vd;vu_ud vu_vd];