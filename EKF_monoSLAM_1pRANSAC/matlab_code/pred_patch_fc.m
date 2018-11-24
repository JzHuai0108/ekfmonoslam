%-----------------------------------------------------------------------
% 1-point RANSAC EKF SLAM from a monocular sequence
%-----------------------------------------------------------------------

% Copyright (C) 2010 Javier Civera and J. M. M. Montiel
% Universidad de Zaragoza, Zaragoza, Spain.

% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation. Read http://www.gnu.org/copyleft/gpl.html for details

% If you use this code for academic work, please reference:
%   Javier Civera, Oscar G. Grasa, Andrew J. Davison, J. M. M. Montiel,
%   1-Point RANSAC for EKF Filtering: Application to Real-Time Structure from Motion and Visual Odometry,
%   to appear in Journal of Field Robotics, October 2010.

%-----------------------------------------------------------------------
% Authors:  Javier Civera -- jcivera@unizar.es 
%           J. M. M. Montiel -- josemari@unizar.es

% Robotics, Perception and Real Time Group
% Aragón Institute of Engineering Research (I3A)
% Universidad de Zaragoza, 50018, Zaragoza, Spain
% Date   :  May 2010
%-----------------------------------------------------------------------

function patch_pred=pred_patch_fc(cam,patch_p_f_ini,R_Wk,r_Wk,XYZ_w)
 
uv_p_pred = patch_p_f_ini.h;

halfW_pred=patch_p_f_ini.half_patch_size_when_matching;

if((uv_p_pred(1)>halfW_pred) && (uv_p_pred(1)<cam.nCols-halfW_pred)&&...
        (uv_p_pred(2)>halfW_pred) && (uv_p_pred(2)<cam.nRows-halfW_pred))  
    
    uv_p_f=patch_p_f_ini.uv_when_initialized;
    R_Wk_p_f=patch_p_f_ini.R_wc_when_initialized;
    r_Wk_p_f = patch_p_f_ini.r_wc_when_initialized;
    
    H_Wk_p_f = [R_Wk_p_f zeros(3,1);zeros(1,3) 1]*[eye(3) r_Wk_p_f; zeros(1,3) 1];
%     H_Wk_p_f = [R_Wk_p_f r_Wk_p_f;zeros(1,3) 1];
    
    H_Wk = [R_Wk zeros(3,1);zeros(1,3) 1]*[eye(3) r_Wk; zeros(1,3) 1];
%     H_Wk = [R_Wk r_Wk;zeros(1,3) 1];
    
    H_kpf_k = inv(H_Wk_p_f)*H_Wk;
    patch_p_f=patch_p_f_ini.patch_when_initialized;
    halfW_fea=patch_p_f_ini.half_patch_size_when_initialized;
    d = cam.dx;
    f = cam.f;
    
    n1 = -[-(uv_p_f(1)-cam.Cx) -(uv_p_f(2)-cam.Cy) f/d]';
    n2 = -[-(uv_p_pred(1)-cam.Cx) -(uv_p_pred(2)-cam.Cy) f/d]';
    n2 = H_kpf_k*[n2;1];
    n2 = n2/n2(4);
    n2 = n2(1:3);
    n1 = n1/norm(n1);
    n2 = n2/norm(n2);
    n = n1+n2;
    n = n/norm(n);
    
    XYZ_kpf = inv(H_Wk_p_f)*[XYZ_w;1];
    XYZ_kpf = XYZ_kpf/XYZ_kpf(4);
    d = -n'*XYZ_kpf(1:3);
    
    uv_p_pred_patch=rotate_with_dist_fc_c2c1(cam,uv_p_f,H_kpf_k(1:3,1:3),H_kpf_k(1:3,4),n,d);
    
    [u_pred,v_pred]=meshgrid(uv_p_pred_patch(1)-halfW_pred:uv_p_pred_patch(1)+halfW_pred,uv_p_pred_patch(2)-halfW_pred:uv_p_pred_patch(2)+halfW_pred);
    uv_pred=[reshape(u_pred,(halfW_pred*2+1)^2,1),reshape(v_pred,(halfW_pred*2+1)^2,1)];
    uv_pred_imak_dist=rotate_with_dist_fc_c1c2(cam,uv_pred,H_kpf_k(1:3,1:3),H_kpf_k(1:3,4),n,d);
    uv_pred_imak_dist(:,1)=uv_pred_imak_dist(:,1)-(uv_p_f(1)-halfW_fea-1);% + 0.5*ones(size(uv_pred_imak_dist,1),1);
    uv_pred_imak_dist(:,2)=uv_pred_imak_dist(:,2)-(uv_p_f(2)-halfW_fea-1);% + 0.5*ones(size(uv_pred_imak_dist,1),1);
    u_pred_imak_dist=reshape(uv_pred_imak_dist(:,1),halfW_pred*2+1,halfW_pred*2+1);
    v_pred_imak_dist=reshape(uv_pred_imak_dist(:,2),halfW_pred*2+1,halfW_pred*2+1);
    [u_fea,v_fea]=meshgrid(1:size(patch_p_f,1),1:size(patch_p_f,2));
    
    patch_pred=interp2(u_fea,v_fea,double(patch_p_f),u_pred_imak_dist,v_pred_imak_dist);
    
else
    patch_pred=zeros(2*halfW_pred+1,2*halfW_pred+1);
end