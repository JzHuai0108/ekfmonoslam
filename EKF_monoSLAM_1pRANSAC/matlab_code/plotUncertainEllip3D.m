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

function plotUncertainEllip3D( C, nu, chi2, color, alpha )

[ X, Y, Z ] = sphere;

nPoints = size(X,2)*size(X,1);

X = reshape( X, 1, nPoints )*sqrt(chi2);
Y = reshape( Y, 1, nPoints )*sqrt(chi2);
Z = reshape( Z, 1, nPoints )*sqrt(chi2);

K = chol( C )';

XYZprima = K*[ X; Y; Z ] + [ ones(1,nPoints)*nu(1); ones(1,nPoints)*nu(2); ones(1,nPoints)*nu(3)];

X = reshape( XYZprima(1,:), sqrt(nPoints), sqrt(nPoints) );
Y = reshape( XYZprima(2,:), sqrt(nPoints), sqrt(nPoints) );
Z = reshape( XYZprima(3,:), sqrt(nPoints), sqrt(nPoints) );

if alpha == 1
    mesh( X, Y, Z, 'FaceColor', color, 'EdgeColor', color, 'LineWidth', 1.5 );
else
    mesh( X, Y, Z, 'FaceColor', color , 'FaceAlpha', 0.2, 'EdgeColor', color, 'EdgeAlpha', 0.4);
end
hold on;