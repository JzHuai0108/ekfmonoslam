%input is defined in "b", ref is defined in "n" (e.g inp=[g*sin(45) g*cos(45) 0], ref=[0 0 1])
%Cbn_est=(I-S(e))Cbn, e=E[inp_err]

function [Cnb E]=align_wander(inp, ref)

%Compute the dcm matrix
Cnb=align_vector(inp,ref);

if (nargout==2) %Compute the error matrix 
    E=zeros(3);
    ep=0.001;
    
    for i=1:3
        %ith perturbation
        inp_err=inp;
        inp_err(i)=inp_err(i)+ep;
        Cnb_err=align_vector(inp_err,ref);
        err=Cnb_err*Cnb'-eye(3);
        E(:,i)=[-err(2,3);err(1,3);-err(1,2)]/ep;
    end
    
end






function Cnb=align_vector(inp, ref)
%normalize
inp=inp/norm(inp);
ref=ref/norm(ref);

%angles
dotp=dot(ref,inp);
crosv=cross(ref,inp);
crosp=norm(crosv);


if (crosp==0) %inp and ref is already aligned
    if (dotp>=0) %no rotation is necessary assert(dotp==1)
        Cnb=eye(3);
    else %assert(dotp==-1)
        %Find a vector that is ortogonal to ref
        if (ref(1)^2>ref(2)^2)
            vr_a=[-ref(3);0;ref(1)];
        else
            vr_a=[0;-ref(3);ref(2)];
        end
        vr_a=vr_a/norm(vr_a)*pi;
        
        %rotate around computed orthogonal vector
        Cnb=rot2dcm_v000(vr_a);
    end
else
    Cnb=eye(3)+skew(crosv)+((1-dotp)/crosp^2)*skew(crosv)*skew(crosv);
end

return;
