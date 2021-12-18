function data=readbin(fname,nrow)

fid=fopen(fname,'rb');
if (fid==-1)
    disp('file not found');
    data=[];
    return;
end

data=fread(fid,[nrow inf],'double');

vr_a=fread(fid,1,'double');
if (~isempty(vr_a))
    disp('readbin: nrow may be wrong!');
end

fclose(fid);