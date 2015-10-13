%Computes the curvature matrix and the gravity.
%Type 1: Non-singular implementation for wander azimuth mechanization
%Type 2: Singular implementation for local geodetic frame mechanization

%Inputs
%Type=1 or 2
%Cen3: 3rd column Earth to nav transformation matrix.
%       [cos(ALPHA)cos(LAT) -sin(ALPHA)cos(LAT) -sin(LAT)]'
%height: Height above elipsoid
%Vn: Velocity in navigation frame

function [Fc, wen_n, wie_n , grav]=geoparam(type, Cen3, height, Vn)

SM_AXIS=6378137;
E_SQR=0.00669437999014;

%Local radii
sr_e=sqrt(1.0-E_SQR*Cen3(3)^2);
Re_h=(SM_AXIS/sr_e)+height;
Rn_h=(SM_AXIS*(1-E_SQR)/(sr_e^3))+height;

if (type==1)
    %Wander azimuth frame mechanization. Note Fc(3,:) will be zero in this implementation
    %Curvature matrix in navigation frame
    sr_a=SM_AXIS*E_SQR/(sr_e^3)/(Re_h)/(Rn_h);
    Fc=zeros(3);
    Fc(1,1)=sr_a*(Cen3(2)^2)+(1/(Re_h));
    Fc(1,2)=-sr_a*Cen3(1)*Cen3(2);
    Fc(2,1)=Fc(1,2);
    Fc(2,2)=sr_a*(Cen3(1)^2)+(1/(Re_h));
    
elseif (type==2)
    %local geodetic frame implementation. Note that, this implementation
    %explicitly uses 1/cosPHI. Therefore, it is singular at the pole.
    
    %Curvature matrix
    Fc=zeros(3);
    Fc(1,1)=1/Re_h;
    Fc(2,2)=1/Rn_h;
    tanLAT=-Cen3(3)/Cen3(1);
    Fc(3,1)=-tanLAT/Re_h;
end


if (nargout>1)
    wen_n=Fc*[Vn(2);-Vn(1);0];  %Note: In wander-azimuth mechanization, wen_n(3)=0
    wie_n=Cen3*7292115e-11;
end

if (nargout==4)
    %Plump bob gravity
    %(Just as a different flavour, I am going to use here a very old gravity model instead of WGS84)
    %See:Heiskanen and Moritz (1967)
    sL2=Cen3(3)^2;
    sL4=Cen3(3)^4;

    a1=9.7803267715;
    a2=0.0052790414;
    a3=0.0000232718;
    a4=-0.0000030876910891;
    a5=0.0000000043977311;
    a6=0.0000000000007211;
    grav=a1*(1+a2*sL2+a3*sL4)+(a4+a5*sL2)*height+a6*height^2;
    
%     %WGS84 approximation for plumb bob gravity.
%     NORMAL_GRV=9.7803253359;
%     GRV_CONS=0.00193185265241;
%     FLATTENING=0.00335281066475;
%     M_FAKTOR=0.00344978650684;
% 
%     sLat=-Cen3(3);
%     g1=NORMAL_GRV*(1+GRV_CONS*sLat*sLat)/sqrt(1.0-E_SQR*sLat^2);
%     grav=g1*(1.0-(2.0/SM_AXIS)*(1.0+FLATTENING+M_FAKTOR-2.0*FLATTENING*sLat*sLat)*height+3.0*height*height/SM_AXIS/SM_AXIS);
end