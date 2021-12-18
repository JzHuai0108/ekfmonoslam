%%Backward computation with BF smoother based on predicted forward
%%solutions. The previous script (smoother_BF_bkwd.m) was based on filtered
%%solution. 

fclose all;
clear;

iniLlh=[39.88864*pi/180;32.78002*pi/180;1112];    %(just for readeble output. not necessary for anything)
DIR='/home/yigiter/Desktop/Dropbox/INS/SAGE/mfiles/data/output/';
nst=21;

%Open the input files to read
fupdres=fopen([DIR 'updresult.bin'],'r');
fpreres=fopen([DIR 'preresult.bin'],'r');
fnavres=fopen([DIR 'prenavresult.bin'],'r');

%Output files
fsmores=fopen([DIR 'smoresult.bin'],'w');

sz_nav=23;
sz_pre=(nst^2+nst*(nst+1)+1);
sz_upd=nst^2+nst*(nst+1)/2+1+nst;

dx_accum=zeros(nst,1);
dP_accum=zeros(nst);

%read the first record to determine when to stop
f_nav=fread(fnavres, sz_nav,'double');
bof=f_nav(1);

%Go to end of file
fseek(fupdres, 0, 'eof');
fseek(fpreres, 0, 'eof');
fseek(fnavres, 0, 'eof');

%%read the last results
fseek(fnavres, -sz_nav*8, 'cof');
fseek(fpreres, -sz_pre*8, 'cof');
fseek(fupdres, -sz_upd*8, 'cof');
f_nav=fread(fnavres,sz_nav,'double');
f_pre=fread(fpreres,sz_pre,'double');
f_upd=fread(fupdres,sz_upd,'double');
proc_time=f_nav(1);
upd_time=f_upd(1);
disp_time=proc_time;

while (1) %I will break inline. (Does octave support exception handling?)
    Llh=f_nav(2:4);
    Vn=f_nav(5:7);
    qbn=f_nav(8:11);
    imuerrors=f_nav(12:end);
    sr_a=(nst*(nst+1)/2);
    P_pre=mat2vec_v000(f_pre(2:1+sr_a));
    STM_pre=reshape(f_pre(2+sr_a:1+sr_a+nst^2),nst,nst);
    Qd=mat2vec_v000(f_pre(2+sr_a+nst^2:end));
    
    %Disp time
    if (disp_time-proc_time(1))>60
        disp(['Process Time:' num2str(proc_time/60)]);
        disp_time=proc_time;
    end
    
    
    if (upd_time==proc_time) %we have observations. Let us correct the predicted soltuion with this observations
        dx_inp=f_upd(2:nst+1);
        dP_inp=mat2vec_v000(f_upd(2+nst:1+nst+(nst*(nst+1)/2)));
        STM_upd=reshape(f_upd(2+nst+(nst*(nst+1)/2):end),nst,nst);
        
        dx_accum=STM_upd'*dx_accum+dx_inp;
        dP_accum=STM_upd'*dP_accum*STM_upd+dP_inp;
                
        %Next update values
        fseek(fupdres, -sz_upd*8*2, 'cof');
        f_upd=fread(fupdres,sz_upd,'double');
        if (isempty(f_upd))
            upd_time=-1;
        else
            upd_time=f_upd(1);
        end
    end
    
    %%Compute the smoothed solution and its variance
    dx=P_pre*dx_accum;
    Llh=Llh-dx(1:3);
    Vn=Vn-dx(4:6);
    qet=rvec2quat_v000(dx(7:9));
    qbn=quatmult_v000(qet, qbn);
    imuerrors=imuerrors+dx(10:end);

    Psmo=P_pre-P_pre*dP_accum*P_pre';
    
    %%Estimated noise
    w_est=Qd*dx_accum;  %note that this noise belongs to (k-1)th epoch when we are at the (k)th epoch
    
    %%Write smoothed solution to the file
    vr_a=posdiff_v000(Llh(1:2)', Llh(3), iniLlh);
    vr_c=quat2dcm_v000(qbn);
    fwrite(fsmores,[proc_time;vr_a';Vn;dcm2euler_v000(vr_c)*180/pi;diag(Psmo)],'double');

    
    %Propagate the backward estimates to next output epoch
    dx_accum=STM_pre'*dx_accum;
    dP_accum=STM_pre'*dP_accum*STM_pre;
    
    
    %Next predicted values
    fseek(fnavres, -sz_nav*8*2, 'cof');
    fseek(fpreres, -sz_pre*8*2, 'cof');
    f_nav=fread(fnavres,sz_nav,'double');
    f_pre=fread(fpreres,sz_pre,'double');
    
    if (f_nav(1)==bof)
        break;
    else
        proc_time=f_nav(1);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Compare results
filtered=readbin_v000([DIR, 'navresult.bin'], 31);
smoothed=readbin_v000([DIR, 'smoresult.bin'], 31);

figure;
plot(filtered(1,:), filtered(2,:));
hold on
plot(smoothed(1,:), smoothed(2,:),'r');

figure;
plot(filtered(2,:), filtered(3,:));
hold on;
plot(smoothed(2,:), smoothed(3,:),'r');

figure;
plot(filtered(1,:), filtered(10+1,:));
hold on
plot(smoothed(1,:), smoothed(10+1,:),'r');
