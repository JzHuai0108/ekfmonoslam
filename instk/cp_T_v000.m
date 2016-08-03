function [T1 Mls]=cp_T(M, R, mtd)
if (size(M,1)>size(M,2))
	[u,s,v]=svd(M);
	T1=u(:,size(M,2)+1:size(M,1))';

	%ls matrix
	if (mtd==0)
		R_inv=inv(R);
		Mls=inv(M'*R_inv*M)*M'*R_inv;
	elseif (mtd==1)
		T2=u(:,1:size(M,2))';   
		mx_a=inv(T1*R*T1');
		mx_b=inv(T2*M);
		Mls=mx_b*(T2-T2*R*T1'*mx_a*T1);
	end
elseif (size(M,1)==size(M,2))
	T1=[]; %No redundancy
	R_inv=inv(R);
	Mls=inv(M'*R_inv*M)*M'*R_inv;
else
	disp('Error cp_T:M is not full col rank');
end
	