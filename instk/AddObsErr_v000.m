function [rstat]=AddObsErr_v000(finname, nrow, foutname, act_row, errdef, rstat)

%randomize
if (rstat(1)~=0)
    randn('state',rstat)
else
    rstat=randn('state');
end

finp=fopen(finname,'rb');
inp=fread(finp,[nrow,inf],'double');
fclose(finp);

for in=1:size(act_row,1)
    fout=fopen(foutname(1,:),'wb');
    act_in=act_row(in,:);
    nd=size(act_in,2);
    if (~isempty(errdef))
        sR=errdef(in).sR;
        for (in1=1:size(inp,2))
            out=inp(act_in,in1)+sR*randn(nd,1);
            fwrite(fout,[inp(1,in1);out],'double');
        end
    else
        out=inp([1 act_in],:);
        fwrite(fout,out,'double');
    end
    fclose(fout);
end
