%Computes the position difference in meters between two coordinates defined
%with Llh.
%Inputs:
%ref=[Lat, lon, height]' of reference.
%height=[h1;h2;...;hn] where n is the number of coordinate
%hpos must be oneof the following: (each position must be on a new row)
%hpos=3rd row of Cen=[-cos(Lat)cos(lon) -cos(Lat)sin(lon) -sin(Lat)];
%or
%hpos=[Lat lon]

%Output: Metric difference between corrdinates defined in the nav frame of
%reference.

function metricdif=posdiff(hpos, height, ref)
E_SQR=0.00669437999014;
SM_AXIS=6378137;

%Convert ref into ecef and compute Cen(ref)
sL=sin(ref(1));
cL=cos(ref(1));
sl=sin(ref(2));
cl=cos(ref(2));
h=ref(3);
    
Re=SM_AXIS/(sqrt(1.0-E_SQR*sL*sL));
refECEF=[(Re+h)*cL*cl (Re+h)*cL*sl (Re*(1-E_SQR)+h)*sL];
CenRef=[-sL*cl -sL*sl cL;-sl cl 0;-cL*cl -cL*sl -sL];

%%%%Convert hpos to ecef
%Be sure that hpos and height is in correct order
if (size(height,2)>size(height,1))
    height=height';
    hpos=hpos';
end

%Compute the normal to the elipsoid in ECEF frame
if size(hpos,2)==2 %hpos is defined as lat/lon. Convert it to 3rd row of Cen
    normal=[-cos(hpos(:,1)).*cos(hpos(:,2)) -cos(hpos(:,1)).*sin(hpos(:,2)) -sin(hpos(:,1))];
elseif (size(hpos,2)==3) %hpos is defined 
    normal=hpos;
end
Re=SM_AXIS./(1.0-E_SQR*normal(:,3).^2).^0.5;

%Compute ecef coordinates
posECEF=-[(Re+height).*normal(:,1) (Re+height).*normal(:,2) (Re*(1-E_SQR)+height).*normal(:,3)];

%Difference
metricdif=(posECEF-ones(size(posECEF,1),1)*refECEF)*CenRef';
