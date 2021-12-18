%etyp==0 => seperate process
%etyp==1 => matrix
function [mod_par sqrt_Q]=exp2markov(exp_par, etyp)
if (etyp==0)    
    nexp=size(exp_par,1);
    mod_par=zeros(nexp,2);
    for i=1:nexp
        a=exp_par(i,1);
        b=(exp_par(i,2)*(1-a^2))^0.5;
        mod_par(i,:)=[a b];
    end
elseif(etyp==1)
    sr_a=size(exp_par,1);
    if (sr_a==1)
        nexp=1;
        Q=exp_par(1,2);
    elseif (sr_a==3)
        nexp=2;
        Q=[exp_par(1,2) exp_par(3,2);exp_par(3,2) exp_par(2,2)];
    elseif (sr_a==6)
        nexp=3;
        Q=[exp_par(1,2) exp_par(4,2) exp_par(6,2);exp_par(4,2) exp_par(2,2) exp_par(5,2);exp_par(6,2) exp_par(5,2) exp_par(3,2)];
    else
        disp('err in exp2markov');
        return;
    end
    
    vr_a=eig(Q);
    a=exp_par(1,1);
    if (~isempty(find(vr_a<0,1,'first')))
        disp('corr error: not psd');
    elseif (~isempty(find(vr_a==0,1,'first')))
        disp('psd');
        sqrt_Q=diag(Q).^0.5;
    else
        sqrt_Q=chol(Q*(1-a^2),'lower');		%%forgot to multiply with the 1-a^2. ahhhhhhhh!!!!
    end
    
    %model for diagonal elements
    mod_par=zeros(nexp,2);
    for i=1:nexp
        a=exp_par(i,1);
        b=(exp_par(i,2)*(1-a^2))^0.5;
        mod_par(i,:)=[a b];
    end
    
end


