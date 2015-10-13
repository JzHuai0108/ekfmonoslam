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

function f = ekf_filter( varargin )

if nargin == 0
    
   f.type = [];
    
   f.x_k_k = [];
   f.p_k_k = [];
   
   f.std_a = [];
   f.std_alpha = [];
   f.sigmaCAM = [];
   
   f.x_k_km1 = [];
   f.p_k_km1 = [];
   f.predicted_measurements = [];
   f.H_predicted = [];
   f.R_predicted = [];
   f.S_predicted = [];
   f.S_matching = [];
   f.z = [];
   f.h = [];
   f.H_matching = [];
   f.measurements = [];
   f.R_matching = [];
   f.x_k_k_mixing_estimate = [];
   f.p_k_k_mixing_covariance = [];
   
   f = class(f,'ekf_filter');
   
elseif isa(varargin{1},'ekf_filter')
    
    f = varargin{1};
    
else
    
   f.type = varargin{6}; 
    
   f.x_k_k = varargin{1};
   f.p_k_k = varargin{2};
   
   f.std_a = varargin{3};
   f.std_alpha = varargin{4};
   f.sigmaCAM = varargin{5};
   
   f.x_k_km1 = [];
   f.p_k_km1 = [];
   f.predicted_measurements = [];
   f.H_predicted = [];
   f.R_predicted = [];
   f.S_predicted = [];
   f.S_matching = [];
   f.z = [];
   f.h = [];
   f.H_matching = [];
   f.measurements = [];
   f.R_matching = [];
   f.x_k_k_mixing_estimate = [];
   f.p_k_k_mixing_covariance = [];
   
%    f = class(f,'ekf_filter');
   
end