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

function [hypothesis_support, positions_li_inliers_id, positions_li_inliers_euc] = compute_hypothesis_support_fast( xi, cam, state_vector_pattern, z_id, z_euc, threshold )

hypothesis_support = 0;

if ~isempty(z_id)
    
    n_id = size(z_id,2);
    
    ri = xi(logical(state_vector_pattern(:,1)));
    anglesi = xi(logical(state_vector_pattern(:,2)));
    rhoi = xi(logical(state_vector_pattern(:,3)));
    
    ri = reshape(ri,3,n_id);
    anglesi = reshape(anglesi,2,n_id);
    mi = m(anglesi);
    
    rwc = xi(1:3);
    rwc = repmat(rwc,1,n_id);
    
    rotwc = q2r(xi(4:7));
    rotcw = rotwc';
    
    ri_minus_rwc = ri - rwc;
    xi_minus_xwc_by_rhoi = ri_minus_rwc(1,:)*diag(rhoi);
    yi_minus_xwc_by_rhoi = ri_minus_rwc(2,:)*diag(rhoi);
    zi_minus_xwc_by_rhoi = ri_minus_rwc(3,:)*diag(rhoi);
    ri_minus_rwc_by_rhoi = [xi_minus_xwc_by_rhoi; yi_minus_xwc_by_rhoi; zi_minus_xwc_by_rhoi];
    
    hc = rotcw*(ri_minus_rwc_by_rhoi + mi);
    
    h_norm = [hc(1,:)./hc(3,:);hc(2,:)./hc(3,:)];
    
    u0 = cam.Cx;
    v0 = cam.Cy;
    f  = cam.f;
    ku = 1/cam.dx;
    kv = 1/cam.dy;
    h_image = f*ku*h_norm + [u0*ones(1,n_id);v0*ones(1,n_id)];
    
    h_distorted = distort_fm( h_image , cam );
    
    nu = z_id - h_distorted;
    residuals = sqrt(nu(1,:).^2+nu(2,:).^2);
    positions_li_inliers_id = residuals<threshold;
    hypothesis_support = hypothesis_support + sum(positions_li_inliers_id);
    
else
    
    positions_li_inliers_id = [];
    
end

if ~isempty(z_euc)
    
    n_euc = size(z_euc,2);
    
    xyz = xi(logical(state_vector_pattern(:,4)));
    xyz = reshape(xyz,3,n_euc);
    
    rwc = xi(1:3);
    rwc = repmat(rwc,1,n_euc);
    
    rotwc = q2r(xi(4:7));
    rotcw = rotwc';
    
    xyz_minus_rwc = xyz - rwc;
    
    hc = rotcw*xyz_minus_rwc;
    
    h_norm = [hc(1,:)./hc(3,:);hc(2,:)./hc(3,:)];
    
    u0 = cam.Cx;
    v0 = cam.Cy;
    f  = cam.f;
    ku = 1/cam.dx;
    kv = 1/cam.dy;
    h_image = f*ku*h_norm + [u0*ones(1,n_euc);v0*ones(1,n_euc)];
    
    h_distorted = distort_fm( h_image , cam );
    
    nu = z_euc - h_distorted;
    residuals = sqrt(nu(1,:).^2+nu(2,:).^2);
    positions_li_inliers_euc = residuals<threshold;
    hypothesis_support = hypothesis_support + sum(positions_li_inliers_euc);
    
else
    
    positions_li_inliers_euc = [];
    
end