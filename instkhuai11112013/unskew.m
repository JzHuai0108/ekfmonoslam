% typ=0, unskew(mx) mx is the common sense skew symmetric matrix of vector
% [1,2,3], i.e., 
% [ 0, -3, 2;
%  3, 0, -1;
% -2, 1, 0]
function [res err]=unskew(mx, typ)
if nargin==1||typ==0
    res=[mx(3,2);mx(1,3);mx(2,1)];
elseif typ==1
    res=-[mx(2,3);mx(3,1);mx(1,2)];
else
    res=([mx(3,2);mx(1,3);mx(2,1)]-[mx(2,3);mx(3,1);mx(1,2)])/2;
end
err=([mx(3,2);mx(1,3);mx(2,1)]+[mx(2,3);mx(3,1);mx(1,2)])./([mx(3,2);mx(1,3);mx(2,1)]-[mx(2,3);mx(3,1);mx(1,2)])*100;