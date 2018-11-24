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

function [Q,G] = func_Q( Xv,u,Pn,delta_t, type )

if strcmp(type,'constant_position_and_orientation_location_noise')
    
    omegaOld=Xv(11:13);
    qOld=Xv(4:7);

    G=sparse(zeros(13,6));

    G(1:3,1:3)=eye(3)*delta_t;
    G(4:7,4:6)=dq_by_deuler(tr2rpy(q2tr(qOld)));
    
else

    omegaOld=Xv(11:13);
    qOld=Xv(4:7);
    qwt=v2q(omegaOld*delta_t);

    G=sparse(zeros(13,6));

    G(8:10,1:3)=eye(3);
    G(11:13,4:6)=eye(3);
    G(1:3,1:3)=eye(3)*delta_t;
    G(4:7,4:6)=dq3_by_dq1(qOld)*dqomegadt_by_domega(omegaOld,delta_t);

end

Q=G*Pn*G';