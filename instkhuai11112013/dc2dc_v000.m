%mode==0 -> c2d
%mode==1 -> d2c
%mode==2 -> d2d
function [At, Qt, Rt]=dc2dc(A,Q,R,dt,mode,dtn)
%%%achtung achtung: code contains recursion. be careful when changing the
%%%version numbers!!!!!!!!!!!!!!!
nst=size(A,1);
Rt=R;
if (mode==0)    %cont 2 disc
    M = dt*[-A,Q;zeros(nst),A'];
    N = expm(M);
    At = N(nst+1:2*nst,nst+1:2*nst)';
    Qt = At*N(1:nst,nst+1:2*nst);
    if (~isempty(R))
        Rt=R/dt;
    end
elseif (mode==1)    %disc 2 cont
    At=logm(A)/dt;      %achtung: eig(A) must not contain "0"
    mx_a=expm(At'*dt);
    mx_b=expm(-At*dt);
    mx_c=mx_b*Q;
    
    M=logm([mx_b mx_c;zeros(nst) mx_a]);
    Qt=M(1:nst,nst+1:2*nst)/dt;
    
    if (~isempty(R))
        Rt=R*dt;
    end
elseif (mode==2) %d2d in discrete time
    n=round(dtn/dt);
    At=A^n;
    Qt=Q;
    for in=2:n
        Qt=A*Qt*A'+Q;
    end
    Rt=R/n;
elseif (mode==3) %d2d via continuous time
    [Ac, Qc, Rc]=dc2dc_v000(A,Q,R,dt,1);
    [At, Qt, Rt]=dc2dc_v000(Ac,Qc,Rc,dtn,0);
end
