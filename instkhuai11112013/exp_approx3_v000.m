function [exp_apx]=exp_approx3(dboun, exp_par,show_res)

%convert i
nind=find(exp_par(:,2)<=0);
exp_par(nind,2)=-1*exp_par(nind,2);

dind1=dboun(1,1):dboun(1,2);
nd1=dboun(1,2)-dboun(1,1)+1;
dind2=dboun(2,1):dboun(2,2);
nd2=dboun(2,2)-dboun(2,1)+1;
dind3=dboun(3,1):dboun(3,2);
nd3=dboun(3,2)-dboun(3,1)+1;
dind4=dboun(4,1):dboun(4,2);
nd4=dboun(4,2)-dboun(4,1)+1;
nd=nd1+nd2+nd3+nd4;

H=zeros(nd,5);
Y=zeros(nd,1);

H(1:nd1,1)=dind1;
H(1:nd1,3)=1;
H((nd1+1):(nd1+nd2),2)=dind2;
H((nd1+1):(nd1+nd2),4)=1;
H((nd1+nd2+1):(nd1+nd2+nd3),1)=dind3;
H((nd1+nd2+1):(nd1+nd2+nd3),5)=1;
H((nd1+nd2+nd3+1):nd,2)=dind4;
H((nd1+nd2+nd3+1):nd,5)=1;

Rinv1=1./(linspace(dboun(1,1),dboun(1,2),length(dind1)));
Rinv2=1./(linspace(dboun(2,1),dboun(2,2),length(dind2)));
Rinv3=1./(linspace(dboun(3,1),dboun(3,2)^2,length(dind3)));
Rinv4=1./(linspace(dboun(4,1),dboun(4,2)^2,length(dind4)));
HR(1:nd1,1)=dind1.*Rinv1;
HR(1:nd1,3)=Rinv1;
HR((nd1+1):(nd1+nd2),2)=dind2.*Rinv2;
HR((nd1+1):(nd1+nd2),4)=Rinv2;
HR((nd1+nd2+1):(nd1+nd2+nd3),1)=dind3.*Rinv3;
HR((nd1+nd2+1):(nd1+nd2+nd3),5)=Rinv3;
HR((nd1+nd2+nd3+1):nd,2)=dind4.*Rinv4;
HR((nd1+nd2+nd3+1):nd,5)=Rinv4;

Y(1:nd1)=log(exp_par(1,2)*exp_par(1,1).^dind1);
Y((nd1+1):(nd1+nd2))=log(exp_par(2,2)*exp_par(2,1).^dind2);
Y((nd1+nd2+1):(nd1+nd2+nd3))=log(exp_par(3,2)*exp_par(3,1).^dind3);
Y((nd1+nd2+nd3+1):nd)=log(exp_par(4,2)*exp_par(4,1).^dind4);

res=inv(HR'*H)*HR'*Y;
exp_apx=[exp(res(1)) exp(res(3));exp(res(2)) exp(res(4));exp(res(1)) exp(res(5));exp(res(2)) exp(res(5))];
exp_apx(nind,2)=-1*exp_apx(nind,2);

if (show_res==1)
    exp_par(nind,2)=-1*exp_par(nind,2);
    figure;
    for i=1:4
        ind=0:dboun(i,2);
        plot(ind,exp_par(i,1).^ind*exp_par(i,2))
        hold on
        plot(ind,exp_apx(i,1).^ind*exp_apx(i,2),'r')
    end
    grid;
end