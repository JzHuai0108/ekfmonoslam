%%Backward computations of RTS smoother
%%Method 2: In this implementation, instead of using a seperate accumulator
%%for smoother corrections I will use smoothed states itself.
%%Compare this implementation with method 1. As you can see they are almost
%%the exactly the same. This implementaion is generally regarded as standard RTS
%%solution. However, as I noted previouly I alwasy prefer BF smoothers (and
%%hence the other RTS implementation) than this.

%%% VERY IMPORTANT NOTE %%%
%%You must note that PROC_A (the procedure to compute dx_inp and PROC_B
%%(the procedure to correct updates states with dx_smo) must be exactly
%%inverse of each other. Otherwise this implementation diverges. Any states
%%other than attitudes are simply added/substracted. Therefore,
%%this inverse operation condition is not an issue for them. However, for
%%attitude states this is a huge problem. As an example uncomment the code
%%parts in the PROC_A and PROC_B, and rerun this code. Although (under
%%small angle assumption) these uncommented parts can also be assumed
%%inverse of each other, you will see that the smoother result diverge.
%%This is one of the reason why you should not use this implementation as
%%RTS smoother. (Now, you may understand why BF is my favorite smoother)
%%(Also you must note that FS is a stable matrix (compute its eigenvalues at each cycle and see yourself)
%%Despite its stability the results still diverge. I find this fact very
%%irritating. This point deserves much deeper analysis, yet, I do not have any more patience to
%%deal with this code any further. The explanation of this divergence problem can be
%%a subject of a fantastic paper.)

fclose all;
clear;

iniLlh=[39.88864*pi/180;32.78002*pi/180;1112];    %(just for readable output. not necessary for anything)
DIR='/home/yigiter/Desktop/Dropbox/INS/SAGE/mfiles/data/output/';   %Output Directory of the forward path where the smoother buffer is saved
nst=21;
nhlf=nst*(nst+1)/2;

%Open the input files to read
fsmobuf=fopen([DIR 'smobuffer.bin'],'r');
sz_upd=23+nhlf;
sz_pre=23+nhlf+nst^2;

%Output files
fsmores=fopen([DIR 'smoresult.bin'],'w');

%read the first record to determine when to stop
smobuf=fread(fsmobuf, sz_pre+sz_upd,'double');
bof=smobuf(1);
clear smobuf;

%Go to end of file
fseek(fsmobuf, 0, 'eof');

%%read the last results
fseek(fsmobuf, -(sz_pre+sz_upd)*8, 'cof');
vr_a=fread(fsmobuf,sz_pre+sz_upd,'double');
predicted=vr_a(1:sz_pre);
updated=vr_a(sz_pre+1:end);
disp_time=predicted(1);

%Initialize the smoothed solution
Llh_smo=updated(2:4);
Vn_smo=updated(5:7);
qbn_smo=updated(8:11);
imuerrors_smo=updated(12:23);
P_smo=mat2vec_v000(updated(24:end));

% %%DEBUG
% dx_accum=zeros(21,1);

while (1)
    %%Write the smoothed solution
    vr_a=posdiff_v000(Llh_smo(1:2)', Llh_smo(3), iniLlh);
    vr_c=quat2dcm_v000(qbn_smo);
    fwrite(fsmores,[predicted(1);vr_a';Vn_smo;dcm2euler_v000(vr_c)*180/pi;imuerrors_smo;diag(P_smo)],'double');

    %Disp time
    if (disp_time-predicted(1))>60
        disp_time=predicted(1);
        disp(['Process Time:' num2str(disp_time/60)]);
    end
    
    %Compute the difference between the smoothed and predicted solutions (PROC_A)
    dx_inp=zeros(21,1);
    dx_inp(1:3)=predicted(2:4)-Llh_smo; %note that the sign is different from the usuall convention.!!!!
    dx_inp(4:6)=predicted(5:7)-Vn_smo;
    qpu=quatmult_v001(qbn_smo, predicted(8:11),2);
    dx_inp(7:9)=quat2rot_v000(qpu);
%     mx_a=quat2dcm_v000(qbn_smo);
%     mx_b=quat2dcm_v000(predicted(8:11));
%     mx_c=mx_a*mx_b';
%     dx_inp(7:9)=[mx_c(3,2);mx_c(1,3);mx_c(2,1)];
%     %dx_inp(7:9)=dcm2rot_v000(mx_c);
    dx_inp(10:21)=imuerrors_smo-predicted(12:23);

%     %%%DEBUG%%%%
%     %Compute the difference between the smoothed and predicted solutions
%     dx=zeros(21,1);
%     dx(1:3)=predicted(2:4)-updated(2:4); %note that the sign is different from the usuall convention.!!!!
%     dx(4:6)=predicted(5:7)-updated(5:7);
%     mx_a=quat2dcm_v000(updated(8:11));
%     mx_b=quat2dcm_v000(predicted(8:11));
%     mx_c=mx_a*mx_b';
%     dx(7:9)=[mx_c(3,2);mx_c(1,3);mx_c(2,1)];   
%     dx(10:21)=updated(12:23)-predicted(12:23);
%     oriantation states
    STM_pre=reshape(predicted(24+nhlf:end),nst,nst);
    P_pre=mat2vec_v000(predicted(24:23+nhlf));
    %FS1=STM_pre'/P_pre;
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
    dx_smo=FS*dx_inp;
    
    %%%%DEBUG (assert dx_smo==dx_accum)
%     dx_accum=FS*(dx_accum+dx);
%     %dx_smo=dx_accum;
    
    %correct the filtered states with the smoothed estimates (PROC_B)
    Llh_smo=updated(2:4)-dx_smo(1:3);
    Vn_smo=updated(5:7)-dx_smo(4:6);
%     mx_a=eye(3)+skew(dx_smo(7:9));
%     mx_b=mx_a*quat2dcm_v000(updated(8:11));
%     qbn_smo=dcm2quat_v000(mx_b);
    qet=rvec2quat_v000(dx_smo(7:9));
    qbn_smo=quatmult_v000(qet, updated(8:11));
    imuerrors_smo=updated(12:23)+dx_smo(10:end);
    
    if (predicted(1)==bof)
        break;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Compare results
filtered=readbin_v000([DIR, 'navresult.bin'], 31+12);
smoothed=readbin_v000([DIR, 'smoresult.bin'], 31+12);

figure;
plot(filtered(1,:), filtered(4,:));
hold on
plot(smoothed(1,:), smoothed(4,:),'r');

figure;
plot(filtered(1,:), filtered(11,:));
hold on
plot(smoothed(1,:), smoothed(11,:),'r');

figure;
plot(filtered(1,:), filtered(8,:));
hold on
plot(smoothed(1,:), smoothed(8,:),'r');

figure;
plot(filtered(2,:), filtered(3,:));
hold on;
plot(smoothed(2,:), smoothed(3,:),'r');
% 
% figure;
% plot(filtered(1,:), filtered(10+12+8,:));
% hold on
% plot(smoothed(1,:), smoothed(10+12+8,:),'r');
