%if input is a symmetric matrice, converts it into a vector using diagonal
%order of matrix (main diagonal will be the first elements of the vector).
%if the input is a vector, generates a symmetric matrix from the vector.

function out = mat2vec(inp)

[nrow ncol]=size(inp);
if (nrow==ncol) %inp is a(symmetric) matrix
    out=zeros(nrow*(nrow+1)/2, 1);
    m=1;
    n=nrow;
    for i=1:nrow
        out(m:m+n-1)=diag(inp,i-1);
        m=m+n;
        n=n-1;
    end
else %inp is a vector. convert it to matrix
    n=max(nrow,ncol);
    nrow=(-1+sqrt(1+8*n))/2;
    if ((nrow-round(nrow))~=0)
        disp('Error:smat2vec:21');
        out=0;
    else
        out=zeros(nrow);
        m=nrow+1;
        n=nrow-1;
        out=out+diag(inp(1:nrow));
        for i=1:nrow-1
            out=out+diag(inp(m:m+n-1),i)+diag(inp(m:m+n-1),-i);
            m=m+n;
            n=n-1;
        end
    end
end

        