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
function [scrc,crc] = crosscorr(I1,I2,SVD)

% CROSSCORR Computes the normalized cross-correlation between two images
%
% [scrossc crossc] = crosscorr(I1,I2)
%
% input
%   I1,I2: Images
%
% output
%   scrossc: sum of the cross-correlations
%   crossc:  vector of cross-correlations
%
% NOTE: the correlation of 0_nxn and 0_nxn will be zero

%Check input parameters
sz = ones(1,3);
sz(1:ndims(I1)) = size(I1);
if prod(double((size(I2)==size(I1))))==0, %jmmm because prod didn't multipy logical values
  error(' ');
end;

%The same with singular values (invariant to rotations)
if (nargin==3)&&(SVD=='svd'),
  [scrc,crc] = crosscorrsvd(I1,I2);
  return
end;

%Compute the normalized cross-correlation
flag= 1; %See std for details
num = (I1-repmat(mean(mean(I1,1),2),sz(1:2))).*(I2-repmat(mean(mean(I2,1),2),sz(1:2)));
den = repmat(std(reshape(I1,[1,prod(sz(1:2)),sz(3)]),flag,2),sz(1:2))...
    .*repmat(std(reshape(I2,[1,prod(sz(1:2)),sz(3)]),flag,2),sz(1:2));
crc = (den~=0).*num./(den+(den==0));
scrc= reshape(mean(mean(crc,1),2),[1,size(crc,3)]);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Auxiliar functions

function [scrc,crc] = crosscorrsvd(I1,I2)

flag = 1;
crc = [];
for i=1:size(I1,3),
  d1 = svd(I1(:,:,i)); d1 = diag(d1);
  d2 = svd(I2(:,:,i)); d2 = diag(d2);
  num = (d1-mean(d1,1)).*(d2-mean(d2,1));
  den = repmat(std(d1,flag,1).*std(d2,flag,1),size(num));
  crc = cat(3,crc,(den~=0).*num./(den+(den==0)));
end;
scrc = mean(crc,1);

return


