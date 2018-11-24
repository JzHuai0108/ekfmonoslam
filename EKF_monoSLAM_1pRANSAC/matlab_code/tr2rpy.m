%TR2RPY	Convert a homogeneous transform matrix to roll/pitch/yaw angles
%
%	[A B C] = TR2RPY(TR) returns a vector of Euler angles 
%	corresponding to the rotational part of the homogeneous transform TR.
%
%	See also  RPY2TR, TR2EUL

%	Copright (C) Peter Corke 1993
function rpy = tr2rpy(m)
	
	rpy = zeros(1,3);

	if abs(m(1,1)) < eps & abs(m(2,1)) < eps,
		rpy(1) = 0;
		rpy(2) = atan2(-m(3,1), m(1,1));
		rpy(3) = atan2(-m(2,3), m(2,2));
	else,
		rpy(1) = atan2(m(2,1), m(1,1));
		sp = sin(rpy(1));
		cp = cos(rpy(1));
		rpy(2) = atan2(-m(3,1), cp * m(1,1) + sp * m(2,1));
		rpy(3) = atan2(sp * m(1,3) - cp * m(2,3), cp*m(2,2) - sp*m(1,2));
	end
