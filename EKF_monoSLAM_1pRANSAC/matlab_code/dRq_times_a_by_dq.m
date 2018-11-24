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

function dRq_times_a_by_dqRES=dRq_times_a_by_dq(q,aMat)

  dRq_times_a_by_dqRES=zeros(3,4);
  
  TempR = dR_by_dq0(q);
  Temp31 = TempR * aMat;
  dRq_times_a_by_dqRES(1:3,1)=Temp31;

  TempR = dR_by_dqx(q);
  Temp31 = TempR * aMat;
  dRq_times_a_by_dqRES(1:3,2)=Temp31;

  TempR = dR_by_dqy(q);
  Temp31 = TempR * aMat;
  dRq_times_a_by_dqRES(1:3,3)=Temp31;

  TempR = dR_by_dqz(q);
  Temp31 = TempR * aMat;
  dRq_times_a_by_dqRES(1:3,4)=Temp31;
  
 return

 
function dR_by_dq0RES=dR_by_dq0(q)
  q0 = q(1);
  qx = q(2);
  qy = q(3);
  qz = q(4);

  dR_by_dq0RES=[2*q0, -2*qz,  2*qy;
		        2*qz,  2*q0, -2*qx;
		       -2*qy,  2*qx,  2*q0];
 
  return;


function dR_by_dqxRES=dR_by_dqx(q)

  q0 = q(1);
  qx = q(2);
  qy = q(3);
  qz = q(4);

  
  dR_by_dqxRES=[2*qx, 2*qy,   2*qz;
		        2*qy, -2*qx, -2*q0;
		        2*qz, 2*q0,  -2*qx];
 
return



function dR_by_dqyRES=dR_by_dqy(q)
    
  q0 = q(1);
  qx = q(2);
  qy = q(3);
  qz = q(4);

  dR_by_dqyRES=[-2*qy, 2*qx,  2*q0;
		         2*qx, 2*qy,  2*qz;
		        -2*q0, 2*qz, -2*qy];
 
return

function dR_by_dqzRES=dR_by_dqz(q)
  q0 = q(1);
  qx = q(2);
  qy = q(3);
  qz = q(4);


  dR_by_dqzRES=[-2*qz, -2*q0, 2*qx;
		         2*q0, -2*qz, 2*qy;
		         2*qx,  2*qy, 2*qz];
 
  return
