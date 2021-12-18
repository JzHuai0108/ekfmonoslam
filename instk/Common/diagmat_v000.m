%diag_mat(A,B,2)=[B 0 0;0 A 0;0 0 A];
function [res]=diagmat(A,B,rep)

if (isempty(A))
    res=B;
    return;
end

if (isempty(rep))
    rep=1;
end

nrow=size(A,1);
ncol=size(A,2);

if (isempty(B))  
    res=zeros(nrow*rep,ncol*rep);
    for in=1:rep
        res((in-1)*nrow+1:in*nrow,(in-1)*ncol+1:in*ncol)=A;
    end
else
    nrow1=size(B,1);
    ncol1=size(B,2);
    res=zeros(nrow*rep+nrow1,ncol*rep+ncol1);
    res(1:nrow1,1:ncol1)=B;
    for in=1:rep
        res((in-1)*nrow+1+nrow1:in*nrow+nrow1,(in-1)*ncol+1+ncol1:in*ncol+ncol1)=A;
    end
end