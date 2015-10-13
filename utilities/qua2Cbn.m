function Cbn = qua2Cbn(qua, isnorm);
%-------------------------------------------------------
% C_nb = qua2Cbn(qua, isnorm);
% Yudan Yi, May 26, 2005
%-------------------------------------------------------
% qua = > C_nb
q0 = qua(1);
q1 = qua(2);
q2 = qua(3);
q3 = qua(4);
if (nargin<2) isnorm = true; end;
if (isnorm)
	norm_q = sqrt(q0*q0+q1*q1+q2*q2+q3*q3);
	q0 = q0/norm_q;
	q1 = q1/norm_q;
	q2 = q2/norm_q;
	q3 = q3/norm_q;
end
Cbn(1,1) = q0*q0+q1*q1-q2*q2-q3*q3;
Cbn(1,2) = 2*(-q0*q3+q1*q2);
Cbn(1,3) = 2*( q0*q2+q1*q3);
Cbn(2,1) = 2*( q0*q3+q1*q2);
Cbn(2,2) = q0*q0-q1*q1+q2*q2-q3*q3;
Cbn(2,3) = 2*(-q0*q1+q2*q3);
Cbn(3,1) = 2*(-q0*q2+q1*q3);
Cbn(3,2) = 2*( q0*q1+q2*q3);
Cbn(3,3) = q0*q0-q1*q1-q2*q2+q3*q3;
return;
