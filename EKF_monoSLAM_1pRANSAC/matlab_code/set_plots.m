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
% paths
directory_storage_name = 'figures';
if( exist( directory_storage_name, 'dir' ) )
else
    mkdir( directory_storage_name );
end

figure_all = figure;
oldUnits = get( figure_all, 'Units' );
set( figure_all, 'Units', 'normalized' );
set( figure_all, 'Position', [0.01,0.2,0.7,0.4] );
set( figure_all, 'Units', oldUnits );

% models_fig = subplot('position',[0.05 0.1 0.05 0.8]);
im_fig = subplot('position',[0.05 0.1 0.425 0.8]);
near3D_fig = subplot('position',[0.575 0.1 0.375 0.8]);

view(-180,0);