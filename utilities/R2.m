function N = R2(rad)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% R2 rotation matrix
% WI07:GS660 - Geodetic Reference System
%
% YoungJin Lee 
% Graduate Student
% Geodetic Science & Surveying
% The Ohio State University
% 02-28-2007
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = [cos(rad) 0 -sin(rad);0 1 0;sin(rad) 0 cos(rad)];