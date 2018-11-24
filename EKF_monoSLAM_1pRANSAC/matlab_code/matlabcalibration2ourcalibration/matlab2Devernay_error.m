function error=matlab2Devernay_error(X,camera, uv_d,xy_matlab_N)

% J Huai. Jul 2014
%
% Error between the Devernay undistortion model and the
%   Swaminathan model
% Input 
%   X       -Devernay model params
%   uv_d    -Distorted point coordinates
%   xy_matlab_N    -Undistorted coordinates according to Bouguet model in
%   z=1 plane

camera.omega=X(1);

xy_Devernay_N = undistort_and_normalize_Devernay( uv_d, camera );
if(1) 
    cla;
    plot(xy_matlab_N(1,:),xy_matlab_N(2,:),'+r');
    hold on
    plot(xy_Devernay_N(1,:),xy_Devernay_N(2,:),'+g');
    xlabel('red, matlab.    green Devernay with 1 radial distorion parameter')
end
error=[xy_Devernay_N(1,:)-xy_matlab_N(1,:);xy_Devernay_N(2,:)-xy_matlab_N(2,:)];
end

function xyu=undistort_and_normalize_Devernay( uvd, camera )

xd = ( uvd(1,:) - camera.cx )/camera.fx;
yd = ( uvd(2,:) - camera.cy )/camera.fy;

rd = sqrt( xd.*xd + yd.*yd );
ru= tan(rd * camera.omega)/(2*tan(camera.omega/2));
D= ru./rd;
xyu = [xd.*D;yd.*D];
end