function [Mls, AEq, BEq, CEq, DEq, sPEq]=cp_redsys_v000(M, A, B, C, D, sP)
nax=size(M,2);
nout=size(M,1);
if (nout<=nax)
    Mls=M;
    AEq=A;
    BEq=B;
    CEq=C;
    DEq=D;
    sPEq=sP;
else %redundant observation
    R=D*D'*eye(nout);
    [Mls, Mlserr]=lsmat_v001(M,R,0);
    Rr=Mls*R*Mls';
    DEq=Rr^0.5;
    CEq=(Mls*Mls')^0.5*diagmat_v000(C,[],nax);
%     DEq=chol(Mlserr, 'lower');
%     mx_a=chol(Mls*Mls', 'lower');
%     CEq=mx_a*diagmat_v000(C,[],nax);
    
    AEq=diagmat_v000(A,[],nax);
    BEq=diagmat_v000(B,[],nax);
    sPEq=diagmat_v000(sP,[],nax);
end
