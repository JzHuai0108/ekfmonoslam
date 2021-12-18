function xyz = blh2xyz(blh, a, f);
%--------------------------------------------------------------------------
% xyz = blh2xyz(blh, a, f);
% coordinate transformation from blh to xyz
% note: assume blh is row-based data
% Author: Yudan Yi, Feb. 2005
%--------------------------------------------------------------------------
xyz = [];
if (nargin<1) error('error in data'); end;
if (nargin<2) a = 6378137.0; f = 1.0/298.257223563; end;
if (nargin<3)                f = 1.0/298.257223563; end;
if (isempty(a)) a = 6378137.0; end;
if (isempty(f)) f = 1.0/298.257223563; end;
if (isempty(blh)) error('error in data'); end;
[n,m] = size(blh);
if (n~=3 & m~=3) error('error in data'); end;
if (n==3) 
    if (m~=3) blh = blh'; end;
end
[n,m] = size(blh);
for i=1:n
    lat = blh(i,1); lon = blh(i,2); ht = blh(i,3);
    if (abs(lat)>(2.0*pi)&abs(lon)>(2.0*pi))
        lat = lat*pi/180.0;
        lon = lon*pi/180.0;
    end
    e2 = 2*f-f*f;
    W = sqrt(1-e2*sin(lat)*sin(lat));
    N = a/W;
    xyz(i,1) = (N+ht)*cos(lat)*cos(lon);
    xyz(i,2) = (N+ht)*cos(lat)*sin(lon);
    xyz(i,3) = (N*(1-e2)+ht)*sin(lat);
end