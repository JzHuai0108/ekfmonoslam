function z = geoatan2(x, y)
%-------------------------------------------------------
% z = geoatan2(x, y);
% Yudan Yi, May 26, 2005
%-------------------------------------------------------
% return -pi~pi of y/x
z = 0.0;
if (abs(x)<1e-12)  if (y>0) z = pi/2.0; else z =-pi/2.0; end; return; end;
z = atan(abs(y/x));
if (x>0)  
	if (y<0) z =-z; end; return;
else 
	if (y>0) z = pi-z; else z =-pi+z; end; 
end
return;
