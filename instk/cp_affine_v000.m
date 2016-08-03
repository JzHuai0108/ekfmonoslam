%dout=A*din+b;
function [A b]=cp_affine(dout,din)
nvin=size(din,1);
nvout=size(dout,1);
ndat=size(din,2);

data=[din;ones(1,ndat)];
mx_a=(dout*data')/(data*data');

A=mx_a(1:nvout,1:nvin);
b=mx_a(1:nvout,nvin+1);

