%%Inputs:
%gyro   :angle increment type gyro outputs
%acc    :velocity increment type acc outputs
%mint   :number of minor interval to be used
%alg    :desired coning algorithm

%%outputs
%vinc   :velocity increment in the final period (in a period equal to
%mint*msubint)
%vcorr  :sculling correction over the increment. compensated
%output=vinc+vcorr
%ainc   :angle increment in the final period
%acorr  :coning correction. final compansated output=ainc+acorr

%%Note: if you simply want to downsample some imu data use alg=0 and set
%%mint to the desired downsampling factor. (essentialy this also
%%corresponds to ignagni(1990):algA.)
function [vinc, vcorr, ainc, acorr]=sculling_coning(gyro, acc, mint, alg)

%Use roscoe(2001)'s equivalency rule for ignagni's methods
[ainc, acorr]=coning_v000(gyro, mint, alg);
[vinc, v1corr]=coning_v000(acc, mint, alg);
[v2inc, v2corr]=coning_v000(gyro+acc, mint, alg);

%%Sculling correction
sculling=v2corr-v1corr-acorr;

%%rotation correction
rot=0.5*cross(ainc,vinc,1);

vcorr=sculling+rot;

