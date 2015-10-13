function [exp_apx]=exp_approx(dboun, exp_par, show_res)
nexp=size(exp_par,1);

nind=find(exp_par(:,2)<=0);
exp_par(nind,2)=-1*exp_par(nind,2);

H=zeros(nexp+1);
Y=zeros(nexp+1,1);

for i=1:nexp
    dind=(dboun(i,1):dboun(i,2));
    R=1./(linspace(dboun(i,1),dboun(i,2)^3,length(dind)));
    H(1,1)=H(1,1)+(dind.*R)*dind';
    H(1,1+i)=sum(dind.*R);
    H(1+i,1)=H(1,1+i);
    H(i+1,i+1)=sum(R);
    
    data=dind*log(exp_par(i,1))+log(exp_par(i,2));
    Y(1,1)=Y(1,1)+(dind.*R)*data';
    Y(i+1,1)=sum(R.*data);
end

log_sol=inv(H)*Y;
a=exp(log_sol(1));
K=exp(log_sol(2:end));
K(nind)=-1*K(nind);
exp_apx=[ones(nexp,1)*a K];

if (show_res==1)
    exp_par(nind,2)=-1*exp_par(nind,2);
    figure;
    for i=1:nexp
        ind=0:dboun(i,2);
        plot(ind,exp_par(i,1).^ind*exp_par(i,2))
        hold on
        plot(ind,exp_apx(i,1).^ind*exp_apx(i,2),'r')
    end
    grid;
end