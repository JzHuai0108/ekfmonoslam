function Xsol=Bouguet2Tsai()
%
% JMM Montiel. 6 sep 2005

% load matlab_calibration.mat
% fc = fc_frontal;
% kc = kc_frontal;
% cc = cc_frontal;
% alpha_c = alpha_c_frontal;
% cam = initialize_cam;
% Huai{
load H:\relaylatest\toolbox_calib\calib_casios\Calib_Results.mat
cam.k1 =    -kc(1);
cam.k2 =    kc(2);
cam.nRows = 720;
cam.nCols = 1280;
cam.Cx =    cc(1);
cam.Cy =    cc(2);
d=0.0112;
cam.f =     d*fc(1);
cam.dx =    d;
cam.dy =    d;
cam.model = 'two_distortion_parameters';


cam.K =     sparse( [ fc(1)   0     cc(1);
                0  fc(2)    cc(2);
                0    0     1] );

% Huai}

% % generate distorted table
% nCols = 320;
% nRows = 240;
% inc=round(min(nCols,nRows)/20)
% [u_grid,v_grid]=meshgrid(-150:inc:nCols+150,-150:inc:nRows+150);
% [nGrid,mGrid]=size(u_grid);
% undistorted=[reshape(u_grid,1, nGrid*mGrid); reshape(v_grid,1,nGrid*mGrid)];
% distorted = distort_bouguet(undistorted,kc,fc,cc)

addpath TOOLBOX_Calib;


Xini = [cam.k1; cam.k2; cam.f; cam.Cx; cam.Cy]

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
x_distort = [(xy_u(1,:)*fc(1) + cc(1));(xy_u(2,:)*fc(2) + cc(2))];

% Non-linear minimization, minimizes the error in the image grid between the
% Bouguet's calibration model and the one we use in our code. Xsol should
% return the calibration parameters for our code in the following order:
% k1, k2, f, cx, cy. dx and dy are set here to 0.0112
Xsol = lsqnonlin('matlab2tsai_error',...
    Xini,[0;0;0;0;0],[3*cam.k1;3*cam.k2;3*cam.f;nCols;nRows],...
    optimset('LargeScale','on','Display','Iter','TolFun',1e-14),cam,uv_d,x_distort)