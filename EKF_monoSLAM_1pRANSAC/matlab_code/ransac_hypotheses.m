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

function [features_info] = ransac_hypotheses( filter, features_info, cam )

p_at_least_one_spurious_free = 0.99; % default value
% RANSAC threshold should have a low value (less than the standard
% deviation of the filter measurement noise); as high innovation points
% will be later rescued
threshold = get_std_z(filter);

n_hyp = 1000; % initial number of iterations, will be updated
max_hypothesis_support = 0; % will be updated

[state_vector_pattern, z_id, z_euc] = generate_state_vector_pattern( features_info, get_x_k_km1(filter) );

for i = 1:n_hyp
    
    % select random match
    [zi,position,num_IC_matches] = select_random_match(features_info);
    
    % 1-match EKF state update
    x_k_km1 = get_x_k_km1(filter);
    p_k_km1 = get_p_k_km1(filter);
    hi = features_info(position).h';
    Hi = sparse(features_info(position).H);
    S = full(Hi*p_k_km1*Hi' + features_info(position).R);
    K = p_k_km1*Hi'*inv(S);
    xi = x_k_km1 + K*( zi - hi );
    
    % predict features
    % Compute hypothesis support: predict measurements and count matches
    % under a threshold
    % features_info_i = predict_camera_measurements( xi, cam, features_info );
    % hypothesis_support = count_matches_under_a_threshold( features_info_i );
    % Fast version of computing hypothesis support, avoiding loops with
    % high cost in Matlab
    [hypothesis_support, positions_li_inliers_id, positions_li_inliers_euc] = compute_hypothesis_support_fast( xi, cam, state_vector_pattern, z_id, z_euc, threshold );
    
    if hypothesis_support > max_hypothesis_support
        max_hypothesis_support = hypothesis_support;
        features_info = set_as_most_supported_hypothesis( features_info, positions_li_inliers_id, positions_li_inliers_euc );
        epsilon = 1-(hypothesis_support/num_IC_matches);
        n_hyp = ceil(log(1-p_at_least_one_spurious_free)/log(1-(1-epsilon)));
        if n_hyp==0 break; end
    end
    
    if i>n_hyp break; end
    
end