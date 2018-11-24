function error=Tsai2SwaminathanError(X,camera,uv_d,xy_u)

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
f=X(2);
camera.k1=k1;
camera.f=f;
camera.K(1,3)=X(3);
camera.K(2,3)=X(4);
uv_u_SN=undistort(camera,uv_d);
xy_u_SN=[uv_u_SN(:,1)-camera.K(1,3), uv_u_SN(:,2)-camera.K(2,3)]*camera.dx/f;
if(1) 
    cla;
    plot(xy_u_SN(:,1),xy_u_SN(:,2),'+r');
    hold on
    plot(xy_u(:,1),xy_u(:,2),'+g');
    xlabel('red, Swaminathan model.    green Tsai with 1 radial distorion parameter')
end
error=[xy_u(:,1)-xy_u_SN(:,1);xy_u(:,2)-xy_u_SN(:,2)];