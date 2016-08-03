%(a0)x(n+m)+(a1)x(n+m-1)...+(am)x(n)=(b0)w(n+m)+(b1)w(n+m-1)...+(bk)a(n+m-k) 
%note poles must be <1 and k<m otherwise the results will be wrong
function res=arma2cov(B, A)
m=size(A,2)-1;
k=size(B,2)-1;
m_exp=(m:-1:0)';
k_exp=(k:-1:0)';
[r,p]=residue(B,A);

res_coef=zeros(m,1);
for in=1:m
    res_coef(in)=r(in)* (B*((p(in)^-1).^k_exp)) / (A*((p(in)^-1).^m_exp)) / p(in);
end

res=[p res_coef];
