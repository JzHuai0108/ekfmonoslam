function N = R1(rad)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% R1 rotation matrix
% WI07:GS660 - Geodetic Reference System
%
% YoungJin Lee 
% Graduate Student
% Geodetic Science & Surveying
% The Ohio State University
% 02-28-2007
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = [1 0 0;0 cos(rad) sin(rad);0 -sin(rad) cos(rad)];