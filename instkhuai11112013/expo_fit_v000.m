%data=Ka^ind+M;
%data must start from zero (K+M, Ka+M,....)
%data boun is the first and final ind values to be used (can be 0:end)

%ctyp: the method used to compute expoenential
%ctyp==0 => log based least square (all values must be postivie)
%ctyp==1 => normal equation (can be used also for negative values)

function [dmar_par exp_par]=expo_fit(dboun, data, ctyp, show_res, cval)
[ndim ndat]=size(data);
if (ctyp==2)
    nbas=size(cval,2);
else
    nbas=1;
end
a=zeros(ndim,nbas);
M=zeros(ndim,1);
K=zeros(ndim,nbas);

for i=1:ndim
    lb=dboun(i,1);
    ub=dboun(i,2);
    
    if (ctyp==1)
        %%normal equations based computation
        sr_a=data(i,(lb+1):ub)*data(i,(lb+1):ub)';
        sr_b=data(i,(lb+1):ub)*data(i,(lb+2):ub+1)';
        a_apx=sr_b/sr_a;
        dind=lb:ub;
    elseif (ctyp==0)        
        %%apply a log based least square to compute a
        ind=(lb:ub);
        pind=find(data(i,ind+1)>eps*100);
        dind=ind(pind);

        H=[ones(length(dind),1) dind'];
        Y=log(data(i,dind+1))';
        x=inv(H'*H)*H'*Y;
        K_apx0=exp(x(1));
        a_apx=exp(x(2));
    elseif (ctyp==2)    %exponential bases are provided
        a_apx=cval(i,:);
        dind=lb:ub;
    end
    
    %optimal K & M given suboptimal "a"
    H=zeros(length(dind),nbas+1);
    for k=1:nbas
        H(:,k)=a_apx(k).^dind;
    end
    H(:,nbas+1)=1;
    
    apx_sol=inv(H'*H)*H'*data(i,dind+1)';
    K_apx=apx_sol(1:nbas)';
    M_apx=apx_sol(nbas+1);

    a(i,:)=a_apx;
    K(i,:)=K_apx;
    M(i)=M_apx;

    %compare the results
    if (show_res==1)
        figure;
        %compare the results
        ind=0:(length(data(i,:))-1);
        plot(ind,data(i,:));grid;
        hold on
        val_tot=0;
        for k=1:nbas
            val_th=K(i,k)*a(i,k).^ind+M(i);
            val_tot=val_tot+val_th-M(i);
            plot(ind,val_th,'r');
        end
        plot(ind,val_tot+M(i),'r');
    end
end

%exponential param
exp_par=[a K M];

%discrete markov param x'=ax+bu
b=(K.*(ones(ndim,1)-a.^2)).^0.5;
dmar_par=[a b];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %optimal solution given suboptimal "a"
%     vr_a=a_apx.^dind;
%     sr_a=sum(vr_a);
%     sr_b=data(i,dind+1)*vr_a';
%     sr_c=sum(vr_a.^2);
%     sr_d=sum(data(i,dind+1));
%     sr_e=length(dind);
% 
%     sr_f=sr_a^2-sr_e*sr_c;
%     M_apx=(sr_a*sr_b-sr_c*sr_d)/sr_f;
%     K_apx=(sr_a*sr_d-sr_b*sr_e)/sr_f;
    
%     %apply curve fit
%     %fit function
%     myfun2=@(p,k) p(1)+p(2)*p(3).^k;
%     opt_param=lsqcurvefit(myfun2,[M_apx, K_apx, a_apx],ind,data(i,:));
%     a(i)=opt_param(3);
%     K(i)=opt_param(2);
%     M(i)=opt_param(1);