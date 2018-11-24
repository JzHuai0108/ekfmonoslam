function [uv_u]=TsaiImageUndistort(uv_d,Uo,Vo,dx,kappa)
%
% JMM Montiel 6 sep 2005. 
% Image undistortion for a set of points
%
% input
%  uv_d - distorted point image coordinates.
%         a row per point.
%  Uo,V0,dx,kappa, Tsai distortion parameters
%
% output
%  uv_u - undistorted image point coordinates.
%         a row per point


nPoints=size(uv_d,1);
rd_2=dx*dx*((uv_d(:,1)-Uo).^2+(uv_d(:,2)-Vo).^2);
uv_u=[(uv_d(:,1)-Uo).*(1+kappa*rd_2)+Uo, ...
      (uv_d(:,2)-Vo).*(1+kappa*rd_2)+Vo];