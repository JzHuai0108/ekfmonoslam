% input: uv_c1 N x 2 point array in c1 frame, 
% output: uv_c2 N by 2 matrix in C2 image coordinates
% R_c2c1 the rotation from C1 frame to C2 frame, t_c2c1 the translation
% from c1 frame to c2 frame, it is the position of the camera c1 origin in
% c2 frame, n is the average normal direction of the feature point in the 
% C1 frame, d is the distance from the feature to the c1 optical center
% along the minus n direction
% for normalize() input is 2 x N matrix, output is 2 by N matrix,
% project_points2() input 3 by N matrix, output 2 by N matrix
function uv_c2=rotate_with_dist_fc_c2c1(cam, uv_c1, R_c2c1, t_c2c1, n, d)

xn=normalize(uv_c1',cam.fc,cam.cc,cam.kc,cam.alpha);
xn2=((R_c2c1-(t_c2c1*n'/d)))*[xn;ones(1,size(xn,2))];
uv_c2=project_points2(xn2,zeros(3,1),zeros(3,1),cam.fc,cam.cc,cam.kc,cam.alpha); 
uv_c2=uv_c2';