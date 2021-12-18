%%Inputs:
%gyro   :angle increment type gyro outputs
%acc    :velocity increment type acc outputs
%mint   :number of minor interval to be used
%alg    :desired sculling algorithm

%%outputs
%vinc   :velocity increment in the final period (in a period equal to
%mint*msubint)
%vcorr  :sculling correction over the increment. compensated
%output=vinc+vcorr

%Note: the following functions are thereoretically same as
%sculling_coning_v000

function [vinc, corr]=sculling(gyro, acc, mint, alg)

%%Compute the sculling+rotation correction in each minor interval
[mvinc, mainc, mcorr]=sculling_minor_v000(gyro, acc, alg);

len=size(mvinc,2);
outlen=floor(len/mint);
vinc=zeros(3,outlen);
corr=zeros(3,outlen);

ind=1;
for i=mint:mint:len
    tainc=mainc(:,i-mint+1);
    tvinc=mvinc(:,i-mint+1);
    tcorr=mcorr(:,i-mint+1);
    for k=i-mint+2:i
        tcorr=tcorr+mcorr(:,k)+cross(tainc, mvinc(:,k));
        tainc=tainc+mainc(:,k);
        tvinc=tvinc+mvinc(:,k);
    end
    vinc(:,ind)=tvinc;
    corr(:,ind)=tcorr;
    ind=ind+1;
end