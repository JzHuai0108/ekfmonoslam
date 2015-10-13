function C_en = blh2C_en(blh);
%-------------------------------------------------------
% C_en = blh2C_en(blh);
% Yudan Yi, May 26, 2005
%-------------------------------------------------------
lat = blh(1);
lon = blh(2);
ht = blh(3);
C_en(1,1) =-sin(lat)*cos(lon);
C_en(2,1) =-sin(lat)*sin(lon);
C_en(3,1) = cos(lat);
C_en(1,2) =-         sin(lon);
C_en(2,2) =          cos(lon);
C_en(3,2) = 0.0;
C_en(1,3) =-cos(lat)*cos(lon);
C_en(2,3) =-cos(lat)*sin(lon);
C_en(3,3) =-sin(lat);
return;
