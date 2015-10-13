%typ=0->ecef 2 nav || typ=1->nav 2 ecef (nav=Llh in rad)
% from my observation, ecef2geo is almost as accurate as Borkowski and
% Hirvonen iterative scheme, though no so accurate in height, 1e-4 m.
function [out]=ecef2geo(inp,typ)
%Earth Defs.
E_SQR=0.00669437999014;
SM_AXIS=6378137;
FLATTENING=0.00335281066475;

if (typ==0) %ecef2Llh
    X=inp(1);
    Y=inp(2);
    Z=inp(3);
    
    p=sqrt(X^2+Y^2);
    b=SM_AXIS*(1-FLATTENING);
    ep=((SM_AXIS^2)-(b^2))/b^2;
    theta=atan2(Z*SM_AXIS,p*b);
    
    l=atan2(Y,X);
    L=atan2((Z+ep*b*(sin(theta))^3),(p-E_SQR*SM_AXIS*(cos(theta))^3));
    Re=SM_AXIS/sqrt(1-E_SQR*sin(L)^2);
    h=p/cos(L)-Re;
    
    out=[L;l;h];
else %Llh2ecef
    sL=sin(inp(1));
    cL=cos(inp(1));
    sl=sin(inp(2));
    cl=cos(inp(2));
    h=inp(3);
    
    Re=SM_AXIS/(sqrt(1.0-E_SQR*sL*sL));
    X=(Re+h)*cL*cl;
    Y=(Re+h)*cL*sl;
    Z=(Re*(1-E_SQR)+h)*sL;
    out=[X;Y;Z];
end
