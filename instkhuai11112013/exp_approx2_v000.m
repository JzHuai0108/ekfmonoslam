function [mod_res_apx exp_res_apx]=exp_approx2(dboun, exp_par1, exp_par2,show_res)
nexp=size(exp_par1,1);
exp_res_apx=zeros(nexp,3);
mod_res_apx=zeros(nexp,3);

for i=1:nexp
    dind=(dboun(i,1):dboun(i,2));
    
    Y=(exp_par1(i,1).^dind*exp_par1(i,2)-exp_par2(i,1).^dind*exp_par2(i,2))';
    p1=exp_par1(i,1);
    p2=exp_par2(i,1);
    sr_a=p1/(1-p1^2)/(p1-p2)/(1-p1*p2);
    sr_b=p2/(1-p2^2)/(p2-p1)/(1-p1*p2);
    H=[sr_a*p1.^dind'+sr_b*p2.^dind'];
    
    N_est=inv(H'*H)*H'*Y;
    
    exp_res_apx(i,:)=[p1,p2,sqrt(N_est)];
    mod_res_apx(i,:)=[p1+p2,-p1*p2,sqrt(N_est)];
end


if (show_res==1)
    for i=1:nexp
        dind=0:dboun(i,2);
        val_a=exp_par1(i,1).^dind*exp_par1(i,2);
        val_b=exp_par2(i,1).^dind*exp_par2(i,2);
        p1=exp_res_apx(i,1);
        p2=exp_res_apx(i,2);
        N=exp_res_apx(i,3)^2;
        sr_a=N*p1/(1-p1^2)/(p1-p2)/(1-p1*p2);
        sr_b=N*p2/(1-p2^2)/(p2-p1)/(1-p1*p2);
        val_c=sr_a*p1.^dind'+sr_b*p2.^dind';
        figure;
        plot(dind,val_a-val_b);grid;
        hold on
        plot(dind,val_c,'r');
    end
end