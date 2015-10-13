function out=downsample(inp, nsam, ave)
[m,n]=size(inp);
dim=2;
if (m>n)
	dim=1;
end

mx_a=cumsum(double(inp), dim);
if (dim==1)
    out=diff(mx_a(1:nsam:end,:),1,dim)/ave;
else
    out=diff(mx_a(:,1:nsam:end),1,dim)/ave;
end