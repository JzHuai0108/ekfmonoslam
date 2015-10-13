%lvl is the level to compute white noise power
%typ==1 => input is a white noise sequence
%typ==0 => input is rw sequence
function res=cp_wn_pow(input,lvl,typ)

%filter the imput (high pass filtering)
[input_filt, boun]=lvl_filter_v000(input, lvl, 1, 1);

%compute the ratio
nfilt=2^(lvl-1);
if (typ==0) %RW
    sc=(2*nfilt*nfilt+1)/3/nfilt/4;
elseif(typ==1) %WN
    sc=1/nfilt/2;
end

ndat=boun(2)-boun(1)+1;
spow=(input_filt(:,boun(1):boun(2))*input_filt(:,boun(1):boun(2))')/ndat;
wn_pow=spow/sc;

signal_coef=chol(wn_pow);
res=signal_coef;
