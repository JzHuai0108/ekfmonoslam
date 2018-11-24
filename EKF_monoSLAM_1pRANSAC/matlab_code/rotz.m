%ROTZ	Rotation about Z axis
%
%	ROTZ(theta) returns a homogeneous transformation representing a 
%	rotation of theta about the X axis.
%
%	See also ROTX, ROTY, ROTVEC.

% 	Copyright (C) Peter Corke 1990
function r = rotz(t)
	ct = cos(t);
	st = sin(t);
	r =    [ct	-st	0	0
		st	ct	0	0
		0	0	1	0
		0	0	0	1];
