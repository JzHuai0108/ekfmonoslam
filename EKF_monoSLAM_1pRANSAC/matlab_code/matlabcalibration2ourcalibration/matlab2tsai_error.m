function error=matlab2tsai_error(X,camera,uv_d,xy_matlab_N)

% JMM Montiel. Sept 2005
%
% Error between the Tsia undistortion model and the
%   Swaminathan model
% Input 
%   k1      -Swaminathan distortion model
%   camera  -Camera parameters. camera.k1 is irrelevant
%   uv_d    -Distorted point coordinates
%   uv_u    -Undistorted coordinates according to Tsai model

k1=X(1);
k2=X(2);
f=X(3);
cx = X(4);
cy = X(5);

camera.f = f;
camera.Cx = cx;
camera.Cy = cy;
camera.k1 = k1;
camera.k2 = k2;

xy_tsai_N = undistort_and_normalize_fm( uv_d, camera );
if(1) 
    cla;
    plot(xy_matlab_N(1,:),xy_matlab_N(2,:),'+r');
    hold on
    plot(xy_tsai_N(1,:),xy_tsai_N(2,:),'+g');
    xlabel('red, matlab.    green Tsai with 1 radial distorion parameter')
end
error=[xy_tsai_N(1,:)-xy_matlab_N(1,:);xy_tsai_N(2,:)-xy_matlab_N(2,:)];