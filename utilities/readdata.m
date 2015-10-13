function data = readdata(filename, m)
%-------------------------------------------------------------------------------
% data = readdata(filename, m)
% read binary data with m double
% by Dr. Yudan Yi
%-------------------------------------------------------------------------------
data = [];
if (nargin<2) return; end;
if (isempty(m)|m<=0) return; end;
fdat = fopen(filename, 'rb');
if (fdat==-1) return; end;
data = fread(fdat, 'double');
n = size(data,1);
nn = floor(n/m);
data = reshape(data(1:nn*m), m, nn)';
data = sortrows(data, 1);
loc = find(diff(data(:,1))==0);
data(loc, :) = [];
fclose(fdat);
%-------------------------------------------------------------------------------

