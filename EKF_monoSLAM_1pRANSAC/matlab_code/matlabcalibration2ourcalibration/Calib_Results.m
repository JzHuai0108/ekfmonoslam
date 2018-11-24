% Intrinsic and Extrinsic Camera Parameters
%
% This script file can be directly excecuted under Matlab to recover the camera intrinsic and extrinsic parameters.
% IMPORTANT: This file contains neither the structure of the calibration objects nor the image coordinates of the calibration points.
%            All those complementary variables are saved in the complete matlab data file Calib_Results.mat.
% For more information regarding the calibration model visit http://www.vision.caltech.edu/bouguetj/calib_doc/


%-- Focal length:
fc = [ 391.427320917686473 ; 392.028280196873595 ];

%-- Principal point:
cc = [ 342.245081587138714 ; 238.948758664941181 ];

%-- Skew coefficient:
alpha_c = 0.000000000000000;

%-- Distortion coefficients:
kc = [ -0.322710343395157 ; 0.101773957976404 ; -0.000874980840080 ; -0.000135940496402 ; 0.000000000000000 ];

%-- Focal length uncertainty:
fc_error = [ 6.120190672918460 ; 6.094036763132211 ];

%-- Principal point uncertainty:
cc_error = [ 2.973943133826638 ; 3.720699911432898 ];

%-- Skew coefficient uncertainty:
alpha_c_error = 0.000000000000000;

%-- Distortion coefficients uncertainty:
kc_error = [ 0.008122349887786 ; 0.007226565586295 ; 0.001163519629191 ; 0.001046163323157 ; 0.000000000000000 ];

%-- Image size:
nx = 640/2;
ny = 480/2;


%-- Various other variables (may be ignored if you do not use the Matlab Calibration Toolbox):
%-- Those variables are used to control which intrinsic parameters should be optimized

n_ima = 21;						% Number of calibration images
est_fc = [ 1 ; 1 ];					% Estimation indicator of the two focal variables
est_aspect_ratio = 1;				% Estimation indicator of the aspect ratio fc(2)/fc(1)
center_optim = 1;					% Estimation indicator of the principal point
est_alpha = 0;						% Estimation indicator of the skew coefficient
est_dist = [ 1 ; 1 ; 1 ; 1 ; 0 ];	% Estimation indicator of the distortion coefficients


%-- Extrinsic parameters:
%-- The rotation (omc_kk) and the translation (Tc_kk) vectors for every calibration image and their uncertainties

%-- Image #1:
omc_1 = [ 2.235081e+00 ; -2.180885e+00 ; 1.790347e-01 ];
Tc_1  = [ -8.887203e+01 ; 4.470347e+01 ; 1.085683e+03 ];
omc_error_1 = [ 2.506200e-02 ; 2.171023e-02 ; 5.102044e-02 ];
Tc_error_1  = [ 8.323032e+00 ; 1.033251e+01 ; 1.784640e+01 ];

%-- Image #2:
omc_2 = [ 2.261169e+00 ; -2.056628e+00 ; 2.581312e-01 ];
Tc_2  = [ 3.462046e+02 ; 2.067834e+01 ; 9.481134e+02 ];
omc_error_2 = [ 8.665009e-03 ; 1.183034e-02 ; 2.182694e-02 ];
Tc_error_2  = [ 7.389731e+00 ; 9.264077e+00 ; 1.572740e+01 ];

%-- Image #3:
omc_3 = [ 1.644749e+00 ; -2.353315e+00 ; -1.042833e+00 ];
Tc_3  = [ 3.499125e+02 ; -2.155907e+02 ; 7.613752e+02 ];
omc_error_3 = [ 8.855882e-03 ; 1.522325e-02 ; 2.840129e-02 ];
Tc_error_3  = [ 6.229715e+00 ; 7.771306e+00 ; 1.325044e+01 ];

%-- Image #4:
omc_4 = [ 2.057606e+00 ; -1.935199e+00 ; -1.188595e-01 ];
Tc_4  = [ 4.188913e+02 ; -1.437138e+02 ; 7.804221e+02 ];
omc_error_4 = [ 7.162784e-03 ; 1.217134e-02 ; 1.937235e-02 ];
Tc_error_4  = [ 6.308428e+00 ; 7.898064e+00 ; 1.421558e+01 ];

%-- Image #5:
omc_5 = [ 1.956272e+00 ; -1.689284e+00 ; -8.470042e-02 ];
Tc_5  = [ 4.326686e+02 ; 2.666385e+02 ; 8.158327e+02 ];
omc_error_5 = [ 1.860053e-02 ; 2.064811e-02 ; 3.467126e-02 ];
Tc_error_5  = [ 6.983894e+00 ; 8.589073e+00 ; 1.566392e+01 ];

%-- Image #6:
omc_6 = [ 2.014369e+00 ; -1.821516e+00 ; -3.855789e-01 ];
Tc_6  = [ 4.460849e+02 ; 1.278258e+02 ; 7.406669e+02 ];
omc_error_6 = [ 2.117173e-02 ; 2.874970e-02 ; 4.903335e-02 ];
Tc_error_6  = [ 6.513273e+00 ; 7.664232e+00 ; 1.394375e+01 ];

%-- Image #7:
omc_7 = [ 2.029685e+00 ; -1.848621e+00 ; -2.714211e-01 ];
Tc_7  = [ 3.719215e+02 ; -1.156830e+02 ; 7.952150e+02 ];
omc_error_7 = [ 8.618588e-03 ; 1.281868e-02 ; 1.984295e-02 ];
Tc_error_7  = [ 6.378589e+00 ; 7.958956e+00 ; 1.390616e+01 ];

%-- Image #8:
omc_8 = [ 1.231372e+00 ; -2.180602e+00 ; -1.266094e+00 ];
Tc_8  = [ 5.087794e+02 ; -3.421839e+02 ; 7.258855e+02 ];
omc_error_8 = [ 1.111755e-02 ; 1.893945e-02 ; 2.722516e-02 ];
Tc_error_8  = [ 6.528269e+00 ; 8.152515e+00 ; 1.461025e+01 ];

%-- Image #9:
omc_9 = [ -2.115540e+00 ; 2.149814e+00 ; -2.960400e-01 ];
Tc_9  = [ -2.298568e+02 ; -9.077181e+01 ; 5.591121e+02 ];
omc_error_9 = [ 1.053211e-02 ; 9.337858e-03 ; 1.969061e-02 ];
Tc_error_9  = [ 4.636649e+00 ; 5.560313e+00 ; 9.014424e+00 ];

%-- Image #10:
omc_10 = [ -2.002656e+00 ; 2.125271e+00 ; 2.167484e-01 ];
Tc_10  = [ -2.061390e+02 ; -1.107652e+02 ; 5.501064e+02 ];
omc_error_10 = [ 8.970968e-03 ; 1.187485e-02 ; 1.655004e-02 ];
Tc_error_10  = [ 4.421725e+00 ; 5.532352e+00 ; 8.642549e+00 ];

%-- Image #11:
omc_11 = [ -1.995278e+00 ; 1.953823e+00 ; -4.918801e-01 ];
Tc_11  = [ -2.319614e+02 ; -8.106762e+01 ; 6.244913e+02 ];
omc_error_11 = [ 1.234502e-02 ; 1.224508e-02 ; 2.787896e-02 ];
Tc_error_11  = [ 5.349863e+00 ; 6.155556e+00 ; 9.175196e+00 ];

%-- Image #12:
omc_12 = [ -1.959522e+00 ; 1.967779e+00 ; 1.507934e-01 ];
Tc_12  = [ -2.181946e+02 ; -1.745170e+02 ; 5.705101e+02 ];
omc_error_12 = [ 1.359201e-02 ; 1.574184e-02 ; 1.899838e-02 ];
Tc_error_12  = [ 4.578891e+00 ; 5.896079e+00 ; 9.292022e+00 ];

%-- Image #13:
omc_13 = [ -2.038408e+00 ; 2.028516e+00 ; 4.818459e-01 ];
Tc_13  = [ -7.224286e+01 ; -1.151310e+02 ; 5.010870e+02 ];
omc_error_13 = [ 7.686874e-03 ; 9.866858e-03 ; 1.601938e-02 ];
Tc_error_13  = [ 3.903297e+00 ; 4.970922e+00 ; 7.109112e+00 ];

%-- Image #14:
omc_14 = [ -2.195225e+00 ; 2.167444e+00 ; 2.844547e-01 ];
Tc_14  = [ -5.269525e+01 ; -8.909912e+01 ; 3.782363e+02 ];
omc_error_14 = [ 7.354099e-03 ; 8.693979e-03 ; 1.386248e-02 ];
Tc_error_14  = [ 2.947453e+00 ; 3.720571e+00 ; 5.834405e+00 ];

%-- Image #15:
omc_15 = [ 2.211089e+00 ; -2.155123e+00 ; -9.233942e-02 ];
Tc_15  = [ -5.529030e+01 ; -8.649008e+01 ; 3.564566e+02 ];
omc_error_15 = [ 7.951384e-03 ; 7.757206e-03 ; 1.333925e-02 ];
Tc_error_15  = [ 2.797836e+00 ; 3.501649e+00 ; 5.656043e+00 ];

%-- Image #16:
omc_16 = [ 2.256820e+00 ; -2.076571e+00 ; -1.156491e-02 ];
Tc_16  = [ -7.876417e+00 ; 2.774798e+01 ; 3.826953e+02 ];
omc_error_16 = [ 7.115910e-03 ; 5.760119e-03 ; 1.421097e-02 ];
Tc_error_16  = [ 2.917002e+00 ; 3.633016e+00 ; 6.004792e+00 ];

%-- Image #17:
omc_17 = [ 2.165192e+00 ; -2.010654e+00 ; -2.454352e-01 ];
Tc_17  = [ 1.772483e+02 ; 5.687622e+01 ; 4.448630e+02 ];
omc_error_17 = [ 8.323268e-03 ; 9.430919e-03 ; 1.854121e-02 ];
Tc_error_17  = [ 3.494577e+00 ; 4.374034e+00 ; 7.393726e+00 ];

%-- Image #18:
omc_18 = [ 2.096547e+00 ; -2.108930e+00 ; -2.130080e-01 ];
Tc_18  = [ 1.618530e+02 ; -5.786759e+01 ; 4.181046e+02 ];
omc_error_18 = [ 5.860376e-03 ; 8.199379e-03 ; 1.484120e-02 ];
Tc_error_18  = [ 3.252542e+00 ; 4.109848e+00 ; 7.020184e+00 ];

%-- Image #19:
omc_19 = [ 2.311051e+00 ; -2.068002e+00 ; 1.577864e-01 ];
Tc_19  = [ 1.292470e+02 ; -2.925680e+01 ; 4.104170e+02 ];
omc_error_19 = [ 5.266671e-03 ; 7.252421e-03 ; 1.481398e-02 ];
Tc_error_19  = [ 3.117586e+00 ; 3.973329e+00 ; 6.617698e+00 ];

%-- Image #20:
omc_20 = [ 2.397755e+00 ; -1.982549e+00 ; 9.864288e-02 ];
Tc_20  = [ 1.041285e+02 ; 9.170124e+01 ; 4.131072e+02 ];
omc_error_20 = [ 7.376126e-03 ; 6.808112e-03 ; 1.608616e-02 ];
Tc_error_20  = [ 3.163456e+00 ; 3.986278e+00 ; 6.716151e+00 ];

%-- Image #21:
omc_21 = [ 2.358212e+00 ; -1.575553e+00 ; 9.497851e-01 ];
Tc_21  = [ 6.351585e+01 ; 1.946202e+02 ; 5.103947e+02 ];
omc_error_21 = [ 7.467396e-03 ; 8.471958e-03 ; 1.396536e-02 ];
Tc_error_21  = [ 3.982527e+00 ; 4.962941e+00 ; 8.601858e+00 ];


