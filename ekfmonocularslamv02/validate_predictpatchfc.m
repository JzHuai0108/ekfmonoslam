function validate_predictpatchfc()
addpath('H:\relaylatest\instk'); % imu functions
addpath('H:\relaylatest\voicebox\'); % for rotro2qr and rotqr2eu, they are more 
addpath('H:\relaylatest\toolbox_calib');% to use project_point2.m and normalize.m
patch=imread('H:\relaylatest\validatepatchfc.bmp');
patch=double(patch(:,:,1));
figure (1)
imshow(uint8(patch));
poses1=[1,0,0,0,0,0,0]';
qTs2c=poses1;

cam.nCols=10000;
cam.nRows=10000;
cam.fc=[10000;10000];
cam.cc=[5000;5000]; cam.kc=zeros(5,1);
cam.alpha=0;
cam.K=[10000,0,5000; 0,10000,5000;0,0,1];
xync1=[0;0];
rho=1;
% case 1
% poses2=[rotro2qr(R3(45/180*pi)); 0;0;0];
% tempX=R3(45/180*pi)*[0;0;1]/rho;
% tempX=tempX/tempX(3);
% case 2
t=10;
alp=90-atan2(1,t)*180/pi;
poses2=[rotro2qr(R2(-alp/180*pi)); t;0;0];
tempX=R2(-alp/180*pi)*[eye(3), -[t;0;0]]*[0;0;1/rho;1];
tempX=tempX/tempX(3);
uvdc2=cam.K*tempX;
patchc2=pred_patch_fc(cam,xync1,uvdc2(1:2),poses1,poses2, qTs2c, rho, patch);
figure(2)
imshow(uint8(patchc2));
end
function Pcj2cjmk=getRTcj2cjmk(qTs0nsj, qTs0nsjmk, qTs2c)
% the camera's pose c(j-k) and c(j) in s0 frame
qs02cjmk=quatmult_v001(qTs2c(1:4),qTs0nsjmk(1:4), 0);
qs02cj=quatmult_v001(qTs2c(1:4),qTs0nsj(1:4), 0);
% the relation between the frame c(j) and the incomiing
% frame c(j-k)
qcj2cjmk=quatmult_v001(qs02cjmk, qs02cj,2);
Tcj2cjmk=qTs2c(5:7)+quatrot_v000(qs02cjmk,qTs0nsj(5:7)-qTs0nsjmk(5:7), 0)-quatrot_v000(qcj2cjmk, qTs2c(5:7),0);
Pcj2cjmk=[quat2dcm_v000(qcj2cjmk),Tcj2cjmk];
end


% given two camera poses and its appearance in the first camera frame where the
% world coordinate system sits, predict the patch appearance in the second
% camera, the predicted patch. We denote the first camera frame as c1 that
% is the world frame in this function domain, and the second camera pose, c2 
% input: cam, camera intrinsic and distortion parameters, xync1, the 
% normalized image coordinates of the point in c1 frame, 2x1, uvdc2, 2x1, the 
% distorted image coordinates of the point in c2 frame, poses1, qs0 2 s1, 
% T s1 in s0, s0 is a earth fixed sensor frame, poses2, q s0 2 s2 and T s2
% in s0, qTs2c, the sensor frame and the camera frame relation, qs2c and T
% s in c, rho is the inverse depth of the point in c1 frame, patchc1, the
% appearance of the point in the c1 frame
% the homography transform between two point patches are referred to the
% Ch.13 of Hartley and Zisserman multiple view geometry
function patch_pred=pred_patch_fc(cam,xync1,uvdc2,poses1,poses2, qTs2c, rho, patchc1)
half_patch_size_when_initialized=199;
half_patch_size_when_matching=half_patch_size_when_initialized;
halfW_pred=half_patch_size_when_matching;
halfW_fea=half_patch_size_when_initialized;
if((uvdc2(1)>halfW_pred) && (uvdc2(1)<cam.nCols-halfW_pred+1)&&...
        (uvdc2(2)>halfW_pred) && (uvdc2(2)<cam.nRows-halfW_pred+1))
    % note if a minus sign is added before the two 1's in the 
    % following two lines as done in Civera's this function, 
    % This sign change seems to improve the performance a little bit, but
    % it brings about a bug with the line
    % uv_p_pred_patch=rotate_with_dist_fc_c2c1_v001(); 
    % which sometimes gives huge predicted image coordinates and the 
    % following meshgrid cannot guarantee the correct grid size with large numbers
    Pc12c2=getRTcj2cjmk(poses1, poses2, qTs2c);  
    uvdc1=project_points2([xync1;1],zeros(3,1),zeros(3,1),cam.fc,cam.cc,cam.kc,cam.alpha); 
    n1 = [xync1;1];
    n2 = cam.K\[uvdc2;1];
    n2 = Pc12c2(1:3,1:3)'*n2;
    n1 = n1/norm(n1);
    n2 = n2/norm(n2);
    n = n1+n2;
    n = n/norm(n);
    % we put the c1 as the fundamental frame for now 
    d = -n'*[xync1;1]/rho; 
    % the predicted uv image coordinates in c2 frame
    uv_p_pred_patch=rotate_with_dist_fc_c2c1_v001(cam,uvdc1',Pc12c2(1:3,1:3),Pc12c2(1:3,4), n,d);
    % if a bug occurs with the following two lines, it means the geometry
    % of the camera w.r.t the environment is mistaken. Double check it!
    [u_pred,v_pred]=meshgrid(uv_p_pred_patch(1)-halfW_pred:uv_p_pred_patch(1)+halfW_pred,uv_p_pred_patch(2)-halfW_pred:uv_p_pred_patch(2)+halfW_pred);
    uv_pred=[reshape(u_pred,(halfW_pred*2+1)^2,1),reshape(v_pred,(halfW_pred*2+1)^2,1)];
    % the reprojected image coordinates in c1 frame
    uv_pred_imak_dist=rotate_with_dist_fc_c1c2_v001(cam,uv_pred,Pc12c2(1:3,1:3),Pc12c2(1:3,4), n,d);
    % align the predicted patch in c2 to the intensity patch source in c1 frame
    uv_pred_imak_dist(:,1)=uv_pred_imak_dist(:,1)-(uvdc1(1)-halfW_fea-1);% + 0.5*ones(size(uv_pred_imak_dist,1),1);
    uv_pred_imak_dist(:,2)=uv_pred_imak_dist(:,2)-(uvdc1(2)-halfW_fea-1);% + 0.5*ones(size(uv_pred_imak_dist,1),1);
    u_pred_imak_dist=reshape(uv_pred_imak_dist(:,1),halfW_pred*2+1,halfW_pred*2+1);
    v_pred_imak_dist=reshape(uv_pred_imak_dist(:,2),halfW_pred*2+1,halfW_pred*2+1);
    % interpolation to determine the intensity of the patch in c2 according
    % to its c1 corresnpondences
    [u_fea,v_fea]=meshgrid(1:size(patchc1,1),1:size(patchc1,2));
    patch_pred=interp2(u_fea,v_fea,double(patchc1),u_pred_imak_dist,v_pred_imak_dist);    
else
    patch_pred=[];
end
end
