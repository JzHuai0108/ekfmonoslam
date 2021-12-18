%%%C-Type implementation for imu interpolation for iphone data
%oper=new output period
%cdef= col definitions [total, Time , x , y, z]
function [ini_time]=interp_imu(INFILES, OUTFILE, cdef, oper)

% fclose all;
% clear;
% DIRNAME='C:\Documents and Settings\Yigiter\Desktop\data\CalibTest1\';
% INFILES=[[DIRNAME '07-16-00-12-acclog.bin ']; [DIRNAME '07-16-00-12-gyrolog.bin']];
% OUTFILE=[DIRNAME 'imu.bin'];
% cdef=[5 2 3 4 5;5 2 3 4 5];
% oper=2;

%i/o files
nfiles=size(cdef,1);
ndat=size(cdef,2)-1;

infids=zeros(nfiles,1);
for (in=1:nfiles)
    infids(in,1)=fopen(INFILES(in,:),'rb');
end
ofid=fopen(OUTFILE,'wb');


%first data from all files
data=zeros(ndat,nfiles);
for in=1:nfiles
    vr_a=fread(infids(in),[cdef(in,1) 1],'double');
    data(:,in)=vr_a(cdef(in,2:end),1);
end

%Latest data
ct=data(1,1);
for in=2:nfiles
    if ct<data(1,in)
        ct=data(1,in);
    end
end

%Read all files till they past the latest
flag=ones(nfiles,1);
while (sum(flag,1))
    for in=1:nfiles
        if (data(1,in)<ct)
            vr_a=fread(infids(in),[cdef(in,1) 1],'double');
            data(:,in)=vr_a(cdef(in,2:end),1);
        else
            flag(in)=0;
        end
    end
end

ini_time=ct;

%%Write the first data (useless)
oin=0;
vr_a=reshape(data(2:end,:), (ndat-1)*nfiles, 1);
fwrite(ofid, [oin;vr_a],'double');

%Start the main part
ipt=ones(nfiles,1)*ct;
idt=zeros(nfiles,1);
data_inc=zeros(ndat,nfiles);
iflag=zeros(nfiles,1);
fflag=zeros(nfiles,1);


while (~sum(fflag,1))
    %%process all files till they exceed the next sampling point
    for in=1:nfiles
        if (~iflag(in))
            dt=data(1,in)-ipt(in);
            if (idt(in)+dt<=oper)
                data_inc(:,in)=data_inc(:,in)+dt*data(:,in);
                idt(in)=idt(in)+dt;
                
                ipt(in)=data(1,in);
                vr_a=fread(infids(in),[cdef(in,1) 1],'double');
                if (isempty(vr_a))
                    break;  %Why does not feof() work properly! whyyyyyyy
                end
                data(:,in)=vr_a(cdef(in,2:end),1);
            else
                data_inc(:,in)=data_inc(:,in)+(oper-idt(in))*data(:,in);
                ipt(in)=ipt(in)+(oper-idt(in));
                iflag(in)=1;
            end
        end
    end
    
    %Check if one the files reaches the feof
    for in=1:nfiles
        fflag(in)=feof(infids(in));
    end
    
    %%In all files pass the next sampling point, write the results
    if sum(iflag,1)==nfiles
        %%Write data
        oin=oin+1;
        vr_a=reshape(data_inc(2:end,:), (ndat-1)*nfiles, 1)/oper;
        fwrite(ofid, [oin;vr_a],'double');
        
        %Prepare for the next cycle
        idt=zeros(nfiles,1);
        data_inc=zeros(ndat,nfiles);
        iflag=zeros(nfiles,1);
    end
end

%%Close all the files
for in=1:nfiles
    fclose(infids(in));
end
fclose(ofid);

return;