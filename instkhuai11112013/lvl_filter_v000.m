%typ==0 => LP
%typ==1 => HP
%btyp=boundary type (1==only outputs which are not affected by edge effect,
%0==all data)
function [output, boun]=lvl_filter(data, lvl, typ, btyp)
[ndim ndat]=size(data);

if (lvl==0) %no filtering
    output=data;
    boun=[1;length(data)];
    return;
end

nflt=2^(lvl-1); %half of the filter size

if (2*nflt>ndat) %output will be the mean
    if (typ==0)
        dat_mean=mean(data,2);
        output=diag(dat_mean)*ones(ndim,ndat);
        boun=[1;ndat];
    elseif (typ==1)
        output=zeros(ndim,ndat);
        boun=[1;ndat];
    end
    return;
end

heff=ones(1,nflt);
k=0:(2*ndat-1);
W=exp(-sqrt(-1)*pi/(2*ndat)*k);
if (typ==0)
    heff_dct=real(W.*fft(heff,ndat*2,2));
else
    heff_dct=imag(W.*fft(heff,ndat*2,2))*sqrt(-1)*-1;
end
data_dct=real((ones(ndim,1)*W).*fft(data,ndat*2,2));

out_dct=data_dct.*(ones(ndim,1)*heff_dct);
output=4*ifft(out_dct,ndat*2,2,'symmetric')/2/nflt;


%% trim output using the boundary type
if (btyp==0)     %%all data
    output=output(:,1:ndat);
    boun=[nflt+1;ndat-nflt+1];  %index of the valid outputs (time will be ind-1);
elseif (btyp==1)   %%data which are not affected by edge effects
    output=output(:,nflt+1:ndat-nflt+1);
    boun=[1;ndat-2*nflt+1];
end