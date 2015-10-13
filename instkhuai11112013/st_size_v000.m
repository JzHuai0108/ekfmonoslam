function [nst nst_TI nst_TV]=st_size(errdefs)

nsen=length(errdefs);

nst=0;
nst_TV=0;
nst_TI=0;
for in=1:nsen
    nst_TI=nst_TI+size(errdefs(in).A,1);
end


for in=1:nsen
    if (~isempty(errdefs(in).tparam))
        nst_TV=nst_TV+2;
    end
end

nst=nst_TV+nst_TI;