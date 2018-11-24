%ROTX	Rotation about X axis
%
%	ROTX(theta) returns a homogeneous transformation representing a 
%	rotation of theta about the X axis.
%
%	See also ROTY, ROTZ, ROTVEC.

% 	Copyright (C) Peter Corke 1990
function r = rotx(t)
	ct = cos(t);
	st = sin(t);
	r =    [1	0	0	0
		0	ct	-st	0
		0	st	ct	0
		0	0	0	1];
