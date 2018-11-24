function uvu = undistort_fm( uvd, camera )
%
% Undistort image coordinates
% Javier Civera, 16/11/05

% nPoints = size( uvd, 2 );
% uvu = zeros( 2, nPoints );
% for k = 1:nPoints;
%     uvu( :, k ) = undistor_a_point( uvd( :, k ), camera );
% end

Cx = camera.Cx;
Cy = camera.Cy;
k1 = camera.k1;
k2 = camera.k2;
dx = camera.dx / camera.f;
dy = camera.dy / camera.f;

xd = ( uvd(1,:) - Cx )*dx;
yd = ( uvd(2,:) - Cy )*dy;

rd = sqrt( xd.*xd + yd.*yd );

D = 1 + k1*rd.^2 + k2*rd.^4;
xu = xd.*D;
yu = yd.*D;

uvu = [ xu/dx + Cx; yu/dy + Cy ];