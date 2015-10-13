function zi = hi_inverse_depth(hrl, cam)
% zi does not depend on the multiplier rho to the point's coordinates in
% the camera frame, (X^c, Y^c, Z^c)
% hrl= (X^c, Y^c, Z^c)'*rho^c_j

% % Angle limit condition: 
% % If seen 45º from the first time it was seen, do not predict it
% v_corig_p = rho*(features_info.r_wc_when_initialized - t_wc) + mi;
% v_c_p = mi;
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
% is in front of the camera
if(hrl( 3, : )<=0)
    zi=[];
    return;
end

% Is in front of the camera?
if ((atan2( hrl( 1, : ), hrl( 3, : ) )*180/pi < -60) ||...
    (atan2( hrl( 1, : ), hrl( 3, : ) )*180/pi > 60) ||...
    (atan2( hrl( 2, : ), hrl( 3, : ) )*180/pi < -60) ||...
    (atan2( hrl( 2, : ), hrl( 3, : ) )*180/pi > 60))
    zi = [];
    return;
end

% apply distortion Bouguet's approach as introduced in the camera
% calibration toolbox for matlab
%cam.f 2x1 in pixels, fx and fy, c, 2x1 in pixels, principle point coordinates,
%positive, cam.k =[k1,k2,p1,p2,k3] as Bouguet, k3 default to 0, alpha default to 0
uv_d=project_points2(hrl,zeros(3,1),zeros(3,1),cam.fc,cam.cc,cam.kc,cam.alpha); %uv_d 2xN point vector in pixels
% Is visible in the image?
if ( uv_d(1)>0 ) && ( uv_d(1)<cam.nCols ) && ( uv_d(2)>0 ) && ( uv_d(2)<cam.nRows )
    zi = uv_d;
    return;
else
    zi = [];
    return;
end