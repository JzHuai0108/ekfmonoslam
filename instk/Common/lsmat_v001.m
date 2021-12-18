function [LSM LSM_R] =lsmat(H,R, mtd)
if (mtd==0)
    R_inv=inv(R);
    LSM_R=inv(H'*R_inv*H);
    LSM=LSM_R*H'*R_inv;    
elseif (mtd==1)
    [u,s,v]=svd(H);
    T2=u(:,1:size(H,2))';
    T1=u(:,size(H,2)+1:size(H,1))';
    
    mx_a=inv(T1*R*T1');
    mx_b=inv(T2*H);
    LSM=mx_b*(T2-T2*R*T1'*mx_a*T1);
    LSM_R=LSM*R*LSM';
elseif(mtd==2)
    [u,s,v]=svd(R);
    Ht=u'*H;
    Rt=u'*R*u;
    in=find(diag(Rt)<1e-7, 1, 'first');
    if (isempty(in))
        in=size(H,1);
    else
        in=in-1;
    end
    Htr=Ht(1:in,:);
    Rtr_inv=inv(Rt(1:in,1:in));
    LSM_R=inv(Htr'*Rtr_inv*Htr);
    LSM=LSM_R*Htr'*Rtr_inv*u(:,1:in)';    
end
