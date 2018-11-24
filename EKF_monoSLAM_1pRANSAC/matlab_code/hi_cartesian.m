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

function zi = hi_cartesian( yi, t_wc, r_wc, cam, features_info )

% Compute a single measurement
% Javier Civera, 17/11/05

% Points 3D in camera coordinates
r_cw = inv( r_wc );
hrl = r_cw*( yi - t_wc );

% Is in front of the camera?
if ((atan2( hrl( 1, : ), hrl( 3, : ) )*180/pi < -60) ||...
    (atan2( hrl( 1, : ), hrl( 3, : ) )*180/pi > 60) ||...
    (atan2( hrl( 2, : ), hrl( 3, : ) )*180/pi < -60) ||...
    (atan2( hrl( 2, : ), hrl( 3, : ) )*180/pi > 60))
    zi = [];
    return;
end

% % Angle limit condition: 
% % If seen 45º from the first time it was seen, do not predict it
% v_corig_p = yi - features_info.r_wc_when_initialized;
% v_c_p = yi - t_wc;
% alpha = acos(v_corig_p'*v_c_p/(norm(v_corig_p)*norm(v_c_p)));
% if abs(alpha)>pi/4
%     zi = [];
%     return;
% end
% 
% % Scale limit condition:
% % If seen from double or half the scale, do not predict it
% scale = norm(v_corig_p)/norm(v_c_p);
% if (scale>2)||(scale<1/2)
%     zi = [];
%     return;
% end

% Image coordinates
uv_u = hu( hrl, cam );
% Add distortion
uv_d = distort_fm( uv_u , cam );

% Is visible in the image?
if ( uv_d(1)>0 ) && ( uv_d(1)<cam.nCols ) && ( uv_d(2)>0 ) && ( uv_d(2)<cam.nRows )
    zi = uv_d;
    return;
else
    zi = [];
    return;
end