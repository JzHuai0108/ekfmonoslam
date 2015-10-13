function [T1 Mls]=cp_Tparam(M, R)
[u,s,v]=svd(M');
T1=v(:,size(M,2)+1:size(M,1))';

%ls matrix
Mls=lsmat_v000(M,R);