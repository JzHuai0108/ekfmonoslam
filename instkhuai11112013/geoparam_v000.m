%%Given the ned position (L,lon, h) computes the local radious and gravity
function [Rn, Re, gn, sL, cL, WIE_E]=geoparam(pos_n)
%Constants
WIE_E=7292115e-11;    %earth's rotaion rate
SM_AXIS=6378137;
E_SQR=0.00669437999014;
NORMAL_GRV=9.7803253359;
GRV_CONS=0.00193185265241;
FLATTENING=0.00335281066475;
M_FAKTOR=0.00344978650684;


sL=sin(pos_n(1));
cL=cos(pos_n(1));
h=pos_n(3);
Rn=6335439.327292829/(sqrt(1.0-E_SQR*sL*sL)*(1.0-E_SQR*sL*sL));
Re=SM_AXIS/(sqrt(1.0-E_SQR*sL*sL));
g1=NORMAL_GRV*(1+GRV_CONS*sL*sL)/(sqrt(1.0-E_SQR*sL*sL));
gn=g1*(1.0-(2.0/SM_AXIS)*(1.0+FLATTENING+M_FAKTOR-2.0*FLATTENING*sL*sL)*h+3.0*h*h/SM_AXIS/SM_AXIS);
