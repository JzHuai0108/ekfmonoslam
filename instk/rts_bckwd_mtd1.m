%%Backward computations of RTS smoother
%%Method 1: Seperate accumulator for backward corrections. (This is similar
%%to what was implemented in smoother_rts_bkwd except the fact that the
%%dx_input is directly computed form the difference of updated and
%%predicted results.)

fclose all;
clear;

iniLlh=[39.88864*pi/180;32.78002*pi/180;1112];    %(just for readeble output. not necessary for anything)
DIR='/home/yigiter/Desktop/Dropbox/INS/SAGE/mfiles/data/output/';   %Output Directory of the forward path where the smoother buffer is saved 
nst=21;
nhlf=nst*(nst+1)/2;

%Open the input files to read
fsmobuf=fopen([DIR 'smobuffer.bin'],'r');
sz_upd=23+nhlf;
sz_pre=23+nhlf+nst^2;

%Output files
fsmores=fopen([DIR 'smoresult1.bin'],'w');

%read the first record to determine when to stop
smobuf=fread(fsmobuf, sz_pre+sz_upd,'double');
bof=smobuf(1);

%Go to end of file
fseek(fsmobuf, 0, 'eof');

%%read the last results
fseek(fsmobuf, -(sz_pre+sz_upd)*8, 'cof');
vr_a=fread(fsmobuf,sz_pre+sz_upd,'double');
predicted=vr_a(1:sz_pre);
updated=vr_a(sz_pre+1:end);
disp_time=predicted(1);

dx_accum=zeros(21,1);   %%Smoother corrections
P_smo=mat2vec_v000(updated(24:end));
while (1)
    %correct the filtered states with the smoothed estimates to obtain the
    %smoothed values
    Llh_smo=updated(2:4)-dx_accum(1:3);
    Vn_smo=updated(5:7)-dx_accum(4:6);
    qet=rvec2quat_v000(dx_accum(7:9));
    qbn_smo=quatmult_v000(qet, updated(8:11));
    imuerrors_smo=updated(12:23)+dx_accum(10:end);

    %%Write the smoothed solution
    vr_a=posdiff_v000(Llh_smo(1:2)', Llh_smo(3), iniLlh);
    vr_c=quat2dcm_v000(qbn_smo);
    fwrite(fsmores,[predicted(1);vr_a';Vn_smo;dcm2euler_v000(vr_c)*180/pi;imuerrors_smo;diag(P_smo)],'double');

    
    %Disp time
    if (disp_time-predicted(1))>60
        disp_time=predicted(1);
        disp(['Process Time:' num2str(disp_time/60)]);
    end
    
    %Compute the difference between the filtered and predicted solutions to
    %find the dx=P*H*inv(R_innovation)*inno
    dx_inp=zeros(21,1);
    dx_inp(1:3)=predicted(2:4)-updated(2:4); %note that the sign is different from the usuall convention.!!!!
    dx_inp(4:6)=predicted(5:7)-updated(5:7);
    mx_a=quat2dcm_v000(updated(8:11));
    mx_b=quat2dcm_v000(predicted(8:11));
    mx_c=mx_a*mx_b';
    dx_inp(7:9)=[mx_c(3,2);mx_c(1,3);mx_c(2,1)];   
    dx_inp(10:21)=updated(12:23)-predicted(12:23);
    
    STM_pre=reshape(predicted(24+nhlf:end),nst,nst);
    P_pre=mat2vec_v000(predicted(24:23+nhlf));
    % FS1=STM_pre'/P_pre; %Note that after multiplying dx with FS1 we would have STM'*H*inv(R_innovation)*inno just as BF
    R = chol(P_pre);
    FS1=(R\(R'\STM_pre))';
    
    %Read the next data
    status=fseek(fsmobuf, -(sz_pre+sz_upd)*8*2, 'cof');
    vr_a=fread(fsmobuf,sz_pre+sz_upd,'double');
    predicted=vr_a(1:sz_pre);
    updated=vr_a(sz_pre+1:end);
    
    P_upd=mat2vec_v000(updated(24:end));
    FS=P_upd*FS1;
    
    P_smo=FS*(P_smo-P_pre)*FS'+P_upd;
    dx_accum=FS*(dx_accum+dx_inp);
    
    if (predicted(1)==bof)
        break;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Compare results
filtered=readbin_v000([DIR, 'navresult.bin'], 31+12);
smoothed=readbin_v000([DIR, 'smoresult1.bin'], 31+12);

figure;
plot(filtered(1,:), filtered(10,:));
hold on
plot(smoothed(1,:), smoothed(10,:),'r');

figure;
plot(filtered(2,:), filtered(3,:));
hold on;
plot(smoothed(2,:), smoothed(3,:),'r');

figure;
plot(filtered(1,:), filtered(10+12+7,:));
hold on
plot(smoothed(1,:), smoothed(10+12+7,:),'r');
