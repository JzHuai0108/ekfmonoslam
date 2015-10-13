function [dout din]=cp_sensorio(data, bounind, val_set,C)
nobs=size(bounind,1);
nv=size(data,1);

if (isempty(C)) %not a compulsory input
    C=eye(nv);
end

dout=zeros(nv,nobs);
din=zeros(nv,nobs);

for in=1:nobs
    dout(:,in)=mean(data(:,bounind(in,1):bounind(in,2)),2);
    din(:,in)=round2set(C*dout(:,in),val_set);
end


function rvec=round2set(vec, set)
nset=length(set);
nvec=size(vec,1);
rvec=zeros(nvec,1);
for in1=1:nvec
    dval=abs(vec(in1)-set(1));
    val=set(1);
    for in2=1:nset
        if abs(vec(in1)-set(in2))<dval
            val=set(in2);
            dval=abs(vec(in1)-set(in2));
        elseif abs(vec(in1)+set(in2))<dval
            val=-set(in2);
            dval=abs(vec(in1)+set(in2));
        end
    end
    rvec(in1)=val;
end