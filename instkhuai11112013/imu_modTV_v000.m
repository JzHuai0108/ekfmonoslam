function [A, B, C, sP]=imu_modTV(errdefs, tem_dif, tem, iniP)

nsen=length(errdefs);

%form overall system model
A=[];
B=[];
C=[];
sr_a=0;
for in=1:nsen
    if (~isempty(errdefs(in).tparam))
        A=diagmat_v000([1 0;0 1],A,[]);
        B=diagmat_v000([0 0;0 sqrt(abs(tem_dif))*errdefs(in).tparam(2)],B,[]);
        C=diagmat_v000([tem 1],C,[]);
    else
        if (~isempty(C))
            C=[C;zeros(1, size(C,2))];
        else
            sr_a=sr_a+1;
        end
    end 
end
if (sr_a>0 && ~isempty(C))
    C=[zeros(sr_a,size(C,2));C];
end

if (iniP)
    sP=[];
    for in=1:nsen
        if (~isempty(errdefs(in).tparam))
            sP=diagmat_v000([errdefs(in).tparam(1) 0;0 0],sP,[]);
        end
    end
end