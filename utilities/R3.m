function N = R3(rad)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% R3 rotation matrix
% WI07:GS660 - Geodetic Reference System
%
% YoungJin Lee 
% Graduate Student
% Geodetic Science & Surveying
% The Ohio State University
% 02-28-2007
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = [cos(rad) sin(rad) 0;-sin(rad) cos(rad) 0;0 0 1];