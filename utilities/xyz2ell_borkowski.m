function [lat,lon,h]=xyz2ell_borkowski(X,Y,Z,a,b,e2)
% XYZ2ELL3  Converts cartesian coordinates to ellipsoidal.
%   Uses direct algorithm in K.M. Borkowski, 1989, Borkowski, Kazimierz M. 
% "Accurate algorithms to transform geocentric to geodetic coordinates." 
%Bulletin géodésique 63.1 (1989): 50-56.
%   only applicable to northern latitude and further than 45km from the earth
%   center
% Version: 2013-09-19
% Useage:  [lat,lon,h]=xyz2ell3(X,Y,Z,a,b,e2)
%          [lat,lon,h]=xyz2ell3(X,Y,Z)
% Input:   X \
%          Y  > vectors of cartesian coordinates in CT system (m)
%          Z /
%          a   - ref. ellipsoid major semi-axis (m); default GRS80
%          b   - ref. ellipsoid minor semi-axis (m); default GRS80
%          e2  - ref. ellipsoid eccentricity squared; default GRS80
% Output:  lat - vector of ellipsoidal latitudes (radians)
%          lon - vector of ellipsoidal longitudes (radians)
%          h   - vector of ellipsoidal heights (m)

if nargin ~= 3 & nargin ~= 6
  warning('Incorrect number of input arguments');
  return
end
if nargin == 3
  [a,b,e2]=refell('grs80');
end

lon=atan2(Y,X);
r=sqrt(X.*X+Y.*Y);
E=(b*Z-(a^2-b^2))/(a*r);
F=(b*Z+(a^2-b^2))/(a*r);
P=4*(E*F+1)/3;
Q=2*(E^2-F^2);
D=P^3+Q^2;
v=(D^.5-Q)^(1/3)-(D^.5+Q)^(1/3);
G=((E^2+v)^.5+E)/2;
t=sqrt(G^2+(F-v*G)/(2*G-E))-G;
lat=atan(a*(1-t^2)/(2*b*t));
h=(r-a*t)*cos(lat)+(Z-b)*sin(lat);
