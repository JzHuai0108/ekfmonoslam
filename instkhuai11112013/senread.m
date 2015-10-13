%%If the length of the data vector is different for each sensor type, then
%%use the following:

inp = fopen('sensors.bin','r', 'b');
MAX_SEN=15; %15 is assumed maximum number of sensor types
MAX_DAT=9; %9 is assumed maximum data length
out=zeros(MAX_SEN,1);    
data=zeros(MAX_DAT,1); 
while (~feof(inp))
    typ = fread(inp, 1, 'int32',0);
    if (isempty(typ))   %incosistent foef. no more data
        break;
    end
    
    sysTime = fread(inp, 1, 'int64',0);
    evTime = fread(inp, 1, 'int64',0);
    len = fread(inp, 1, 'int32',0);
    for i=1:len
        data(i)=fread(inp, 1, 'float',0);
    end
    
    %Record the data as ascii
    if (out(typ)==0)
        %Create an output file
        out(typ)=fopen(['sensor' num2str(typ) '.txt'],'w');
    end
    
    fprintf(out(typ), '%ld %ld ', sysTime, evTime);
    fprintf(out(typ), ' %.16f ', data(1:len));
    fprintf(out(typ), '\n');
end

fclose(inp);
for i=1:MAX_SEN
    if (out(i)~=0)
       fclose(out(i));
    end
end


% %%In the current version of the android the length of each sensor data vector is 3
% %%(even for the temperature sensor). We can use this fact to avoid loops.
% %%The code below belongs to Thomas Pingel.
% inp = fopen('sensors.bin','r', 'b');
% len=3;
% MAX_SEN=15;
% 
% fseek(inp,0,'bof');
% typ = fread(inp,inf,'int32',(24 + 4*len - 4),'b');
% nRecords = length(typ);
% 
% fseek(inp, 4, 'bof');
% sysTime = fread(inp,inf,'int64',(24 + 4*len - 8),'b');
% 
% fseek(inp, 12, 'bof');
% evTime = fread(inp,inf,'int64',(24 + 4*len - 8),'b');
% 
% data = zeros(nRecords,len);
% for i=1:len
%    fseek(inp, 24 + (i-1)*4, 'bof');
%    data(:,i) = fread(inp,inf,'float',24 + 4*len - 4,'b');
% end
% 
% %Record the data as ascii
% for i=1:MAX_SEN
%     ind=find(typ==i);
%     if (~isempty(ind))
%         dat=[sysTime(ind) evTime(ind) data(ind,:)];
%         save(['sensor' num2str(i) '.txt'], 'dat', '-ASCII','-double');
%     end
% end