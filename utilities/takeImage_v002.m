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
% Arag�n Institute of Engineering Research (I3A)% Universidad de Zaragoza, 50018, Zaragoza, Spain
% Date   :  May 2010
%-----------------------------------------------------------------------

% if downscale is provided, each edge of im is 1/1<<downscale of that of
% the original image
function im=takeImage(image_file_name_prefix,k, image_file_name_suffix, downscale)
if nargin > 2
    imRGB=imread(sprintf('%s%04d.%s',image_file_name_prefix,k,image_file_name_suffix));
else
    imRGB=imread(sprintf('%s%04d.pgm',image_file_name_prefix,k));
end

if size(imRGB, 3) == 3
    im = rgb2gray(imRGB);
else
    im = imRGB(:, :, 1);
end

if nargin == 4
    downorder = 0;
    while downorder < downscale
        im = impyramid(im, 'reduce');
        downorder = downorder + 1;
    end
end
