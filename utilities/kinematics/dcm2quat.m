function qbn = dcm2quat(Cbn)
% note b and n are arbitrary frames
% Yudan, Yi, 2012
s(5) = Cbn(1,1) + Cbn(2,2) + Cbn(3,3);
s(1) = 1.0 + s(5);
s(2) = 1.0 + 2.0*Cbn(1,1)-s(5);
s(3) = 1.0 + 2.0*Cbn(2,2)-s(5);
s(4) = 1.0 + 2.0*Cbn(3,3)-s(5);
qbn=zeros(4,1);
[cmax, imax] = max(s);
switch (imax)
	case 1
		qbn(1) = 0.5*sqrt(s(1));
		qbn(2) = 0.25*(Cbn(3,2) - Cbn(2,3))/qbn(1);
		qbn(3) = 0.25*(Cbn(1,3) - Cbn(3,1))/qbn(1);
		qbn(4) = 0.25*(Cbn(2,1) - Cbn(1,2))/qbn(1);
	case 2
		qbn(2) = 0.5*sqrt(s(2));
		qbn(3) = 0.25*(Cbn(2,1) + Cbn(1,2))/qbn(2);
		qbn(4) = 0.25*(Cbn(1,3) + Cbn(3,1))/qbn(2);
		qbn(1) = 0.25*(Cbn(3,2) - Cbn(2,3))/qbn(2);
	case 3
		qbn(3) = 0.5*sqrt(s(3));
		qbn(4) = 0.25*(Cbn(3,2) + Cbn(2,3))/qbn(3);
		qbn(1) = 0.25*(Cbn(1,3) - Cbn(3,1))/qbn(3);
		qbn(2) = 0.25*(Cbn(2,1) + Cbn(1,2))/qbn(3);
	case 4
		qbn(4) = 0.5*sqrt(s(4));
		qbn(1) = 0.25*(Cbn(2,1) - Cbn(1,2))/qbn(4);
		qbn(2) = 0.25*(Cbn(1,3) + Cbn(3,1))/qbn(4);
		qbn(3) = 0.25*(Cbn(3,2) + Cbn(2,3))/qbn(4);
	otherwise
end
