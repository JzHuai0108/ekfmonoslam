function Cbn = att2Cbn(att);
%--------------------------------------------------------------------------
% calculate the direct cosine matrix c_nb from b => n-frame given att (h,p,r)
% assume 1 row or col input 
% c_nb = r3(-h)*r2(-p)*r1(-r);
% Author:   Yudan Yi
%           Aug. 2004
%           Feb. 2005
%--------------------------------------------------------------------------
if (nargin<1) error('error in data'); end;
att = att(:);
if (length(att)<3) error('error in data'); end;
%H = att(1); P = att(2); R = att(3);
R = att(1); P = att(2); H = att(3);
Cbn(1,1)	= cos(H)*cos(P);
Cbn(2,1)	= sin(H)*cos(P);
Cbn(3,1)	=-	     sin(P);
Cbn(1,2)	=-sin(H)*cos(R)+cos(H)*sin(P)*sin(R);
Cbn(2,2)	= cos(H)*cos(R)+sin(H)*sin(P)*sin(R);
Cbn(3,2)	=                      cos(P)*sin(R);
Cbn(1,3)	= sin(H)*sin(R)+cos(H)*sin(P)*cos(R);
Cbn(2,3)	=-cos(H)*sin(R)+sin(H)*sin(P)*cos(R);
Cbn(3,3)	=                      cos(P)*cos(R);
%--------------------------------------------------------------------------