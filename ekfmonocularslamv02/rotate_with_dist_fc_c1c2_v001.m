% input: uv_c2 N x 2 point array in c2 frame, 
% output: uv_c1 N by 2 matrix in c1 image coordinates
% R_c2c1 the rotation from c1 frame to c2 frame, t_c2c1 the translation
% from c1 frame to c2 frame, it is the position of the camera c1 origin in
% c2 frame, n is the average normal direction of the feature point in the 
% C1 frame, d is the distance from the feature to the c1 optical center
% along the minus n direction
% for normalize() input is 2 x N matrix, output is 2 by N matrix,
% project_points2() input 3 by N matrix, output 2 by N matrix
function uv_c1=rotate_with_dist_fc_c1c2(cam, uv_c2, R_c2c1, t_c2c1, n, d)

xn2=normalize(uv_c2',cam.fc,cam.cc,cam.kc,cam.alpha);
xn=(R_c2c1-(t_c2c1*n'/d))\[xn2;ones(1,size(xn2,2))];
uv_c1=project_points2(xn,zeros(3,1),zeros(3,1),cam.fc,cam.cc,cam.kc,cam.alpha); 
uv_c1=uv_c1';