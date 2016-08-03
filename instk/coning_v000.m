%data   :increment type sensor output
%mint   :number of minor intervals
%alg    :type of coning algorithms to be used (see coning_minor for the
%available algoritms)

function [inc, corr]=coning(data, mint, alg)

%%Compute the coning corrections and increments for each minor interval
[minc mcorr]=coning_minor_v000(data, alg);

len=size(minc,2);
outlen=floor(len/mint);
inc=zeros(3,outlen);
corr=zeros(3,outlen);

ind=1;
for i=mint:mint:len
    tinc=minc(:,i-mint+1);
    tcorr=mcorr(:,i-mint+1);
    for k=i-mint+2:i
        tcorr=tcorr+mcorr(:,k)+0.5*cross(tinc, minc(:,k));
        tinc=tinc+minc(:,k);
    end
    inc(:,ind)=tinc;
    corr(:,ind)=tcorr;
    ind=ind+1;
end