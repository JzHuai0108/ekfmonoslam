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

function [ x_k_k, p_k_k, K ] = update( x_km1_k, p_km1_k, H, R, z, h )

if size(z,1)>0

    % filter gain
    S = full(H*p_km1_k*H' + R);
    K = p_km1_k*H'*inv(S);

    % updated state and covariance
    x_k_k = x_km1_k + K*( z - h );
    p_k_k = p_km1_k - K*S*K';
    p_k_k = 0.5*p_k_k + 0.5*p_k_k';
    % p_k_k = ( speye(size(p_km1_k,1)) - K*H )*p_km1_k;

    % normalize the quaternion
    Jnorm = normJac( x_k_k( 4:7 ) );
    size_p_k_k = size(p_k_k,1);
    p_k_k = [   p_k_k(1:3,1:3)              p_k_k(1:3,4:7)*Jnorm'               p_k_k(1:3,8:size_p_k_k);
        Jnorm*p_k_k(4:7,1:3)        Jnorm*p_k_k(4:7,4:7)*Jnorm'         Jnorm*p_k_k(4:7,8:size_p_k_k);
        p_k_k(8:size_p_k_k,1:3)     p_k_k(8:size_p_k_k,4:7)*Jnorm'      p_k_k(8:size_p_k_k,8:size_p_k_k)];

    x_k_k( 4:7 ) = x_k_k( 4:7 ) / norm( x_k_k( 4:7 ) );

else

    x_k_k = x_km1_k;
    p_k_k = p_km1_k;
    K = 0;

end