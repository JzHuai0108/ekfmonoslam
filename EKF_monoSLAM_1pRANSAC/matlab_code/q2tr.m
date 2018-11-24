%Q2TR	Convert unit-quaternion to homogeneous transform
%
%	T = q2tr(Q)
%
%	Return the rotational homogeneous transform corresponding to the unit
%	quaternion Q.
%
%	See also TR2Q

%	Copyright (C) 1993 Peter Corke
function t = q2tr(q)
	s = q(1);
	x = q(2);
	y = q(3);
	z = q(4);

	r = [	1-2*(y^2+z^2)	2*(x*y-s*z)	2*(x*z+s*y)
		2*(x*y+s*z)	1-2*(x^2+z^2)	2*(y*z-s*x)
		2*(x*z-s*y)	2*(y*z+s*x)	1-2*(x^2+y^2)	];
	t = eye(4,4);
	t(1:3,1:3) = r;
	t(4,4) = 1;
