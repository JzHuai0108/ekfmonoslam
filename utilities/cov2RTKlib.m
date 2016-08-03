% Function converts covariance matrix of 3d coordinates to the RTKlib
% notation of standard deviations and covariances. The conversion is back
% and forth in one function. It supports ecef-XYZ and NED coordinates.
% Type of calculations:
% drc=0 % from 3x3 covariance matrix to 6x1 RTKlib notation
% drc=1 % from 1x6 or 6x1 RTKlib notation to 3x3 covariance matrix

function outp=cov2RTKlib(inp,drc)

if drc==0
    outp=zeros(6,1);
    outp(1:3)=diag(inp).^.5;
    sg=sign(inp);
    outp(4)=sg(1,2)*(abs(inp(1,2)))^.5;
    outp(5)=sg(2,3)*(abs(inp(2,3)))^.5;
    outp(6)=sg(1,3)*(abs(inp(1,3)))^.5;
else
    outp=diag(inp(1:3).^2);
    inp(inp==0)=2*max(abs(inp(:))); % very rarerly, but sometimes var is 0
    sg=sign(inp(4:6));
    outp(1,2)=sg(1)*inp(4)^2; outp(2,1)=outp(1,2);
    outp(2,3)=sg(2)*inp(5)^2; outp(3,2)=outp(2,3);
    outp(1,3)=sg(3)*inp(6)^2; outp(3,1)=outp(1,3);
end