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

function uvd = distort_fm( uv, camera )
%
% Distort image coordinates.
%  The function deals with two models:
%  1.- Real-Time 3D SLAM with Wide-Angle Vision, 
%      Andrew J. Davison, Yolanda Gonzalez Cid and Nobuyuki Kita, IAV 2004.
%  2.- Photomodeler full distortion model.
% input
%    camera -  camera calibration parameters
%    uvd    -  distorted image points in pixels
% output
%    uv     -  undistorted coordinate points 2XN


% nPoints = size( uv, 2 );
% uvd1 = zeros( 2, nPoints );
% for k = 1:nPoints;
%     uvd1( :, k ) = distor_a_point( uv( :, k ), camera );
% end

Cx = camera.Cx;
Cy = camera.Cy;
k1 = camera.k1;
k2 = camera.k2;
dx = camera.dx;
dy = camera.dy;


xu=(uv(1,:)-Cx)*dx;
yu=(uv(2,:)-Cy)*dy;

ru=sqrt(xu.*xu+yu.*yu);
rd=ru./(1+k1*ru.^2+k2*ru.^4);
% rd_km1 = rd;

for k=1:10
    f = rd + k1*rd.^3 + k2*rd.^5 - ru;
    f_p = 1 + 3*k1*rd.^2 + 5*k2*rd.^4;
    rd = rd - f./f_p;
% %     disp(norm(rd-rd_km1))
%     if norm(rd-rd_km1)<100*eps
%         break
%     end
%     rd_km1 = rd;
end

% disp(k);

D = 1 + k1*rd.^2 + k2*rd.^4;
xd = xu./D;
yd = yu./D;

uvd = [ xd/dx+Cx; yd/dy+Cy ];