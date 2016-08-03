%computes the correlation (matrix) sequence given a stable ss description
function [corrmat]=ss2corr_v000(A, Q, C, P0, dt, crind)

nout=size(C,1);
nout=nout*(nout+1)/2;
nst=size(A,1);
corrmat=zeros(nout,length(crind));

%compute the initial covariance
if (isempty(P0))
    %check whether system is stable or not
    eigenval=eig(A);
    if (~isempty(find(eigenval>=1,1)))
        disp('A is not stable');
        eigenval
        return;
    else
        %Use lyapunov eq.
        P0=dlyap(A,Q);
    end
end

%convert the discrete time into cont
crcind=crind*dt;
[Ac, Qc]=dc2dc_v000(A, Q, [], dt, 1,[]);

%compute the covariances
for in=1:length(crind)
    sr_a=crcind(in);
    Ad=expm(Ac*abs(sr_a));
    %[Ad, Qd]=dc2dc_v000(Ac, Qc, [], abs(crcind(in)), 0,[]);
    if (sr_a>0)
        mx_a=C*(Ad*P0)*C';
    elseif (sr_a<0)
        mx_a=C*(P0*Ad')*C';
    else
        mx_a=C*P0*C';
    end
    corrmat(:,in)=mat2vec(mx_a);
end


function vec=mat2vec(mat)
nrow=size(mat,1);
nout=nrow*(nrow+1)/2;
vec=zeros(nout,1);
lb=1;
ub=nrow;
for in=1:size(mat,1)
    vec(lb:ub)=diag(mat,in-1);
    lb=ub+1;
    ub=ub+nrow-in;
end

    
    
