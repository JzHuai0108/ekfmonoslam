function Xsol=Bouguet2Devernay()
%
% J. Huai 19 jul 2014
clear all
close all
load H:\relaylatest\toolbox_calib\calib_casios\Calib_Results.mat
addpath H:\relaylatest\EKF_monoSLAM_1pRANSAC\matlab_code\matlabcalibration2ourcalibration\TOOLBOX_calib;
% for the test case of casio 2

cam.nRows = 720;
cam.nCols = 1280;
cam.fx= fc(1);
cam.fy= fc(2);
cam.cx= cc(1);
cam.cy= cc(2);
Xini =0.99; %-kc(1)*5;

% Image size
nCols = cam.nCols;
nRows = cam.nRows;

% Step for the image grid
inc = round(min(nCols,nRows)/30); 
% Image grid
[u_grid,v_grid] = meshgrid(3*inc:inc:nCols-3*inc,3*inc:inc:nRows-3*inc);
[nGrid,mGrid] = size(u_grid);
uv_d = [reshape(u_grid,1, nGrid*mGrid); reshape(v_grid,1,nGrid*mGrid)];

xy_u = normalize(uv_d,fc,cc,kc,alpha_c);

% Non-linear minimization, minimizes the error in the image grid between the
% Bouguet's calibration model and the one we use in our code. Xsol should
% return the calibration parameters for our code in the following order:
% omega, fx, fy, cx, cy. 
Xsol = lsqnonlin('matlab2Devernay_error',...
    Xini,0,1,...
    optimset('LargeScale','on','Display','Iter','TolFun',1e-9),cam, uv_d,xy_u)