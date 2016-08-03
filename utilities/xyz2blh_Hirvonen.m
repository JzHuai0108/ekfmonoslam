function [blh, iter] = xyz2blh_Hirvonen(xyz, a, e2, eps);
%--------------------------------------------------------------------------
% blh = xyz2blh(xyz, a, e2);
% coordinate transformation from xyz to blh
% note: assume xyz is row-based data
% Author: Yudan Yi, Feb. 2005
%--------------------------------------------------------------------------

if (nargin<1) error('error in data'); end;
if (nargin<2) a = 6378137.0; f = 1.0/298.257223563;e2	= (2-f)*f; end;
if (nargin<3)                f = 1.0/298.257223563;e2	= (2-f)*f; end;
if(nargin==3) eps=1e-12; end
if (isempty(a)) a = 6378137.0; end;
if (isempty(e2)) f = 1.0/298.257223563; e2	= (2-f)*f; end;
if (isempty(xyz)) error('error in data'); end;
[n,m] = size(xyz);
if (n~=3 & m~=3) error('error in data'); end;
if (n==3) 
    if (m~=3) xyz = xyz'; end;
end
[n,m] = size(xyz);
blh =zeros(n,m);
iter =zeros(n,1);
for i=1:n
    x = xyz(i,1); y = xyz(i,2); z = xyz(i,3);
    lon	= geoatan2(x, y);
    dxy	= sqrt(x*x+y*y);
    
    lat0= geoatan2(dxy, z/(1.0-e2));
    lat	= lat0;
    while (true)
        Rw	= sqrt(1.0-e2*sin(lat0)*sin(lat0));
        Rn	= a/Rw;
        lat1= geoatan2(dxy, z+e2*Rn*sin(lat0) );
        iter(i,1)=iter(i,1)+1;
        if (abs(lat1-lat0)<eps) 
            lat = lat1; 
            break; 
        end;
        lat0= lat1;
    end
    Rw	= sqrt(1.0-e2*sin(lat)*sin(lat));
    Rn	= a/Rw;
    if(abs(lat)<pi/180)
    ht	= dxy/cos(lat)- Rn;
    else
        ht	=z/sin(lat)-Rn*(1-e2);
    end
    blh(i,1) = lat; blh(i,2) = lon; blh(i,3) = ht;
end
