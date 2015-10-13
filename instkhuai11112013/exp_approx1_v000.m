function [exp_apx]=exp_approx1(dboun, exp_par1, exp_par2, show_res)
[ndim]=size(exp_par1,1);
a=zeros(ndim,1);
K=zeros(ndim,1);

nind=find(exp_par1(:,2)<=0);
exp_par1(nind,2)=-1*exp_par1(nind,2);

for i=1:ndim
    dind=dboun(i,1):dboun(i,2);
    data=(exp_par1(i,2)*exp_par1(i,1).^dind)-(exp_par2(1,2)*exp_par2(1,1).^dind);
    
    H=[ones(length(dind),1) dind'];
    
    x=inv(H'*H)*H'*log(data)';
    a(i)=exp(x(2));
    K(i)=exp(x(1));
end

K(nind)=-1*K(nind);
exp_apx=[a K];

if (show_res==1)
    exp_par1(nind,2)=-1*exp_par1(nind,2);
    figure;
    for i=1:ndim
        ind=0:dboun(i,2);
        plot(ind,exp_par1(i,1).^ind*exp_par1(i,2)-(exp_par2(1,1).^ind*exp_par2(1,2)))
        hold on
        plot(ind,exp_apx(i,1).^ind*exp_apx(i,2),'r')
    end
    grid;
end