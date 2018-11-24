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

movie = avifile( sprintf( '%s\\results.avi', 'figures' ), 'compression', 'None' );

for i=91:2169
    % if mod(i,2)==0
    fid = fopen(sprintf( '%s\\image%04d.fig', 'figures', i), 'r');
    if (fid~=-1)
        h = openfig(sprintf( '%s\\image%04d.fig', 'figures', i) );
        im = getframe(gcf);
        movie = addframe( movie, im );
        i
        fclose(fid);
        close(h);
    end 
    % end
end

movie = close( movie );