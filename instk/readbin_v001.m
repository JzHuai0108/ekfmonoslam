function data=readbin(fname,nrow,type)
%Open the file
fid=fopen(fname,'rb');
if (fid==-1)
    disp('file not found');
    data=[];
    return;
end

%Determine the type
if (nargin==3)
    if (isempty(type))
        type='double';
    end
else
    type='double';
end

%Read the Data
data=fread(fid,[nrow inf],type);

%Check the last data to see if the nrow is true
sr_a=zeros(1,type);
sr_b=whos('sr_a');
nbyte=sr_b.bytes;

fseek(fid,-nbyte, 'cof');
sr_a=fread(fid,1,type);
if (sr_a~=data(end))
    disp('Achtung: (readbin) "nrow" is probably wrong!');
end

fclose(fid);
