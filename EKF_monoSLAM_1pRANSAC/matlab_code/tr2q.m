%TR2Q	Convert homogeneous transform to a unit-quaternion
%
%	Q = tr2q(T)
%
%	Return a unit quaternion corresponding to the rotational part of the
%	homogeneous transform T.
%
%	See also Q2TR

%	Copyright (C) 1993 Peter Corke
function q = tr2q(t)
	q = zeros(1,4);
	q(1) = sqrt(trace(t))/2;
	kx = t(3,2) - t(2,3);	% Oz - Ay
	ky = t(1,3) - t(3,1);	% Ax - Nz
	kz = t(2,1) - t(1,2);	% Ny - Ox

	if (t(1,1) >= t(2,2)) & (t(1,1) >= t(3,3)) 
		kx1 = t(1,1) - t(2,2) - t(3,3) + 1;	% Nx - Oy - Az + 1
		ky1 = t(2,1) + t(1,2);			% Ny + Ox
		kz1 = t(3,1) + t(1,3);			% Nz + Ax
		add = (kx >= 0);
	elseif (t(2,2) >= t(3,3))
		kx1 = t(2,1) + t(1,2);			% Ny + Ox
		ky1 = t(2,2) - t(1,1) - t(3,3) + 1;	% Oy - Nx - Az + 1
		kz1 = t(3,2) + t(2,3);			% Oz + Ay
		add = (ky >= 0);
	else
		kx1 = t(3,1) + t(1,3);			% Nz + Ax
		ky1 = t(3,2) + t(2,3);			% Oz + Ay
		kz1 = t(3,3) - t(1,1) - t(2,2) + 1;	% Az - Nx - Oy + 1
		add = (kz >= 0);
	end

	if add
		kx = kx + kx1;
		ky = ky + ky1;
		kz = kz + kz1;
	else
		kx = kx - kx1;
		ky = ky - ky1;
		kz = kz - kz1;
	end
	nm = norm([kx ky kz]);
	if nm == 0,
		q = [1 0 0 0];
	else
		s = sqrt(1 - q(1)^2) / nm;
		q(2:4) = s*[kx ky kz];
	end
