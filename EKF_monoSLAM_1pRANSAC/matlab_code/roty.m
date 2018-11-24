%ROTY	Rotation about Y axis
%
%	ROTY(theta) returns a homogeneous transformation representing a 
%	rotation of theta about the Y axis.
%
%	See also ROTX, ROTZ, ROTVEC.

% 	Copyright (C) Peter Corke 1990
function r = roty(t)
	ct = cos(t);
	st = sin(t);
	r =    [ct	0	st	0
		0	1	0	0
		-st	0	ct	0
		0	0	0	1];
