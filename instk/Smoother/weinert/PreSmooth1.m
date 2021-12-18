%% In off-line navigation applications (e.g. post mission trajectory analysis),
%% we may first prefer to process the raw imu data using only no-motion updates.
%% Essentiall, each no motion instant provides an oportunity to self calibrate
%% the imu itself. As smoothers are robust to modelling errors, pre-processing
%% the raw imu data with smoothers for no-motion updates help us to eliminate 
%% unknown (and unmodelled) biases before the imu data are actually processed by the INS.
%% 
%% However, when you use a smoother to generate such a self calibrated imu data,
%% you must also compute the parameters of the Markovian error model parameters
%% of the new (smoothed) data set. (These model parameters are required in navigation
%% filter/smooters).
%% 
%% Computing the smoothing error model parameter can be insanely difficult for
%% the standard smoother implementation such as RTS, BF etc.
%% On the other hand, the intermadiate variables of Weinert's complementary space
%% approach directly corresponds to these smoothing error model parameters. Therefore,
%% weinert's smoothers are perfect for such imu pre-smoothing applications.


%%%%IMU Pre-Smoother
clear;
DIRNAME='C:\Documents and Settings\Yigiter\Desktop\My Dropbox\SRIMU\Smoother\BeniSil';

%For MC 
% rstr=RandStream('mt19937ar','seed',8967);
% RandStream.setDefaultStream(rstr);
rstr=RandStream.getDefaultStream();
reset(rstr,83474);
rstat=rstr.State;


%%Sensor error definitions at 1Hz
%Accelerometer
SenErrDef(1).A=diag([0.9988 1]);
SenErrDef(1).B=diag([3.5e-005 0]);
SenErrDef(1).C=[1 1];
SenErrDef(1).D=diag([0.001]);
SenErrDef(1).sP=diag([0.0007 0.01]);
SenErrDef(1).tparam=[1e-3 5e-5];  %first param=sP of temperature scale factor (RC), second=B of temperature RW

SenErrDef(2).A=diag([0.99885 1]);
SenErrDef(2).B=diag([3.5e-006 0]);
SenErrDef(2).C=[1 1];
SenErrDef(2).D=diag([0.004]);
SenErrDef(2).sP=diag([0.0008 0.001]);
SenErrDef(2).tparam=[4e-3 1e-5];

SenErrDef(3).A=diag([0.9987 1]);
SenErrDef(3).B=diag([3.5e-004 0]);
SenErrDef(3).C=[1 1];
SenErrDef(3).D=diag([0.0001]);
SenErrDef(3).sP=diag([0.00007 0.01]);
SenErrDef(3).tparam=[1e-4 8e-5];

SenErrDef(4).A=diag([0.9989 1]);
SenErrDef(4).B=diag([3.5e-004 0]);
SenErrDef(4).C=[1 1];
SenErrDef(4).D=diag([0.006]);
SenErrDef(4).sP=diag([0.0003 0.06]);
SenErrDef(4).tparam=[1e-4 5e-4];

% SenErrDef(2)=SenErrDef(1);
% SenErrDef(3)=SenErrDef(1);
% SenErrDef(4)=SenErrDef(1);

%Sensor configuration matrix
M=[cos(pi/2) sin(pi/2);
   cos(2*pi/3) sin(2*pi/3);
   cos(pi/4) sin(2*pi/4);
   cos(0) sin(0)];

%temperature model
tem_mod.A=1;
tem_mod.B=0.1;
tem_mod.u=0;

%%%Generate sensor outputs
[rstat]=AddIMUErr_v001([], [DIRNAME '\imu.bin'], [DIRNAME '\imuerr.bin'], 2000, zeros(4,1), SenErrDef, tem_mod, rstat);

%read/create i/o files
[nst nstTI nstTV]=st_size_v000(SenErrDef);
sendat=readbin_v000([DIRNAME '\imu.bin'],size(M,1)+2);
errdat=readbin_v000([DIRNAME '\imuerr.bin'],nst+1);
FSMO_BUF=fopen([DIRNAME '\smobuf.bin'],'wb');
FSMO_RES=fopen([DIRNAME '\smores.bin'],'wb');
FSMO_ERR=fopen([DIRNAME '\smoerr.bin'],'wb');
FSMO_COV=fopen([DIRNAME '\smocov.bin'],'wb');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%just for comparison (not needed)
FFIL_RES=fopen([DIRNAME '\filres.bin'],'wb');
FFIL_ERR=fopen([DIRNAME '\filerr.bin'],'wb');
FFIL_COV=fopen([DIRNAME '\filcov.bin'],'wb');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Sys model for TI part
[Ati, Bti, Cti, Dti, sPti]=imu_modTI_v000(SenErrDef);

%observation model (assumed to be TI for this example)
[T MLS]=cp_T_v000(M,Dti*Dti',0);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Backward iterations %%%%%%%%%%%%%%%%
z=zeros(nst,1);     %backward best estimate=Pinvb*x
Pbinv=zeros(nst);

%Last obs.
y=T*sendat(3:end,end);
tem=sendat(2,end);
[Atv, Btv, Ctv]=imu_modTV_v000(SenErrDef, 0, tem, 0);
Csys=T*[Cti Ctv];
Dsys=T*Dti;
Rinv=inv(Dsys*Dsys');
z=Csys'*Rinv*y;
Pbinv=Csys'*Rinv*Csys;

%Iterate for the rest of the observations
for in=(size(sendat,2))-1:-1:1
    tem=sendat(2,in);
    tem_dif=sendat(2,in+1)-sendat(2,in);
    [Atv, Btv, Ctv]=imu_modTV_v000(SenErrDef, tem_dif, tem, 0);
    
    if (sendat(1,in)>=1000 && sendat(1,in)<=1200)        %Additional IMU observation (in this case, direct observation from stationarity)
        Msys=[M;1 0;0 1];
        Robs=diagmat_v000(zeros(2),Dti*Dti',[]);
        [Tsys MLSsys]=cp_T_v000(Msys,Robs,1);
        yobs=[sendat(3:end,in);zeros(2,1)];
        Cobs=[Cti Ctv;zeros(2,nst)];
    else
        Tsys=T;
        MLSsys=MLS;
        Robs=Dti*Dti';
        yobs=sendat(3:end,in);
        Cobs=[Cti Ctv];
    end
    
    y=Tsys*yobs;    
    Asys=[Ati zeros(nstTI,nstTV);zeros(nstTV,nstTI) Atv];
    Bsys=[Bti zeros(nstTI,nstTV);zeros(nstTV,nstTI) Btv];
    Csys=Tsys*Cobs;
    Rsys=Tsys*Robs*Tsys';

    smo_mat=inv(eye(nst)+Pbinv*Bsys*Bsys');
    Rinv=inv(Rsys);

    %%write to the buffer
    fwrite(FSMO_BUF,[z;smo_mat(:)],'double');

    %compute the new values
    z=Asys*smo_mat*z+Csys'*Rinv*y;
    Pbinv=Asys'*smo_mat*Pbinv*Asys+Csys'*Rinv*Csys;   
end

fclose(FSMO_BUF);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Forward iterations %%%%%%%%%%%%%%%%%%%
%First data (iterations performed on best estimates)
tem=sendat(2,1);
tem_dif=sendat(2,2)-sendat(2,1);
[Atv, Btv, Ctv, sPtv]=imu_modTV_v000(SenErrDef, tem_dif, tem, 1);
Psys=diagmat_v000(sPtv*sPtv',sPti*sPti',[]);
Asys=diagmat_v000(Atv,Ati,[]);
Bsys=diagmat_v000(Btv,Bti,[]);
Qsys=Bsys*Bsys';
Csys=T*[Cti Ctv];
Dsys=T*Dti;
Rsys=Dsys*Dsys';
Rinp=MLS*(Dti*Dti')*MLS';
Cinp=MLS*[Cti Ctv];

%Psmo=inv((inv(Psys)+Pbinv));
Psmo=inv(eye(nst)+Psys*Pbinv)*Psys;
xsmo=Psmo*z;
usmo=MLS*sendat(3:end,1)-Cinp*xsmo;

FSMO_BUF=fopen([DIRNAME '\smobuf.bin'],'rb');
sz_smo=nst*(nst+1);
fseek(FSMO_BUF, -sz_smo*8,'eof');

%err model
vr_a=fread(FSMO_BUF,sz_smo,'double');
z=vr_a(1:nst);
smo_mat=reshape(vr_a(nst+1:end),nst,nst);
Aerr_smo=smo_mat'*Asys; 
Qerr_smo=Qsys*smo_mat;

fwrite(FSMO_RES,[xsmo;usmo],'double');
fwrite(FSMO_ERR,[Aerr_smo(:);mat2vec_v000(Qerr_smo);Cinp(:);mat2vec_v000(Rinp)],'double');
fwrite(FSMO_COV,[mat2vec_v000(Cinp*Psmo*Cinp');mat2vec_v000(Psmo)],'double');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%filtered states (not required; just for comparison)
y=T*sendat(3:end,1);
K=Psys*Csys'*inv(Csys*Psys*Csys'+Rsys);
xfilt=K*y;
ufilt=MLS*sendat(3:end,1)-Cinp*xfilt;
Pfilt=(eye(nst)-K*Csys)*Psys;

Aerr_filt=(Asys-K*Csys*Asys);
Qerr_filt=Qsys+K*Rsys*K';
fwrite(FFIL_RES,[xfilt;ufilt],'double');
fwrite(FFIL_ERR,[Aerr_filt(:);mat2vec_v000(Qerr_filt);Cinp(:);mat2vec_v000(Rinp)],'double');
fwrite(FFIL_COV,[mat2vec_v000(Cinp*Pfilt*Cinp');mat2vec_v000(Pfilt)],'double');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%iterate for the rest
for (in=2:size(sendat,2)-1)   
    %smoothed value
    xsmo=smo_mat'*Asys*xsmo+smo_mat'*Qsys*z;
    Psmo=Aerr_smo*Psmo*Aerr_smo'+Qerr_smo;
    
    %%System model at time "in"
    tem=sendat(2,in);
    tem_dif=sendat(2,in+1)-sendat(2,in);
    [Atv, Btv, Ctv]=imu_modTV_v000(SenErrDef, tem_dif, tem, 0);
    
    if (sendat(1,in)>=1000 && sendat(1,in)<=1200)
        Msys=[M;1 0;0 1];
        Robs=diagmat_v000(zeros(2),Dti*Dti',[]);
        [Tsys MLSsys]=cp_T_v000(Msys,Robs,1);
        yobs=[sendat(3:end,in);zeros(2,1)];
        Cobs=[Cti Ctv;zeros(2,nst)];
    else
        Robs=Dti*Dti';
        Tsys=T;
        MLSsys=MLS;
        yobs=sendat(3:end,in);
        Cobs=[Cti Ctv];
    end
    
    
    Asys=diagmat_v000(Atv,Ati,[]);
    Bsys=diagmat_v000(Btv,Bti,[]);
    Qsys=Bsys*Bsys';
    Csys=Tsys*Cobs;
    Rsys=Tsys*Robs*Tsys';
    Rinp=MLSsys*Robs*MLSsys';
    Cinp=MLSsys*Cobs;
    
    usmo=MLSsys*yobs-Cinp*xsmo;
    
    %error model param
    fseek(FSMO_BUF, -sz_smo*8*2,'cof'); %relocate the file pointer for the next read
    vr_a=fread(FSMO_BUF,sz_smo,'double');
    z=vr_a(1:nst);
    smo_mat=reshape(vr_a(nst+1:end),nst,nst);
    Aerr_smo=smo_mat'*Asys; 
    Qerr_smo=Qsys*smo_mat;
   
    fwrite(FSMO_RES,[xsmo;usmo],'double');
    fwrite(FSMO_ERR,[Aerr_smo(:);mat2vec_v000(Qerr_smo);Cinp(:);mat2vec_v000(Rinp)],'double');
    fwrite(FSMO_COV,[mat2vec_v000(Cinp*Psmo*Cinp');mat2vec_v000(Psmo)],'double');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%filtered solution (just for comparison)
    y=Tsys*yobs;
    K=Pfilt*Csys'*inv(Csys*Pfilt*Csys'+Rsys);
    Pfilt=(eye(nst)-K*Csys)*Pfilt;
    xfilt=xfilt+K*(y-Csys*xfilt);
    ufilt=MLSsys*yobs-Cinp*xfilt;
    
    Aerr_filt=(Asys-K*Csys*Asys);
    Qerr_filt=Qsys+K*Rsys*K';
    fwrite(FFIL_RES,[xfilt;ufilt],'double');
    fwrite(FFIL_ERR,[Aerr_filt(:);mat2vec_v000(Qerr_filt);Cinp(:);mat2vec_v000(Rinp)],'double');
    fwrite(FFIL_COV,[mat2vec_v000(Cinp*Pfilt*Cinp');mat2vec_v000(Pfilt)],'double');
    
    Pfilt=Asys*Pfilt*Asys'+Qsys;
    xfilt=Asys*xfilt;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%last data
tem=sendat(2,end);
tem_dif=0;
[Atv, Btv, Ctv]=imu_modTV_v000(SenErrDef, tem_dif, tem, 0);
Rinp=MLS*(Dti*Dti')*MLS';
Cinp=MLS*[Cti Ctv];

xsmo=smo_mat'*Asys*xsmo+smo_mat'*Qsys*z;
Psmo=Aerr_smo*Psmo*Aerr_smo'+Qerr_smo;
usmo=MLS*sendat(3:end,in)-Cinp*xsmo;
fwrite(FSMO_RES,[xsmo;usmo],'double');
fwrite(FSMO_ERR,[zeros(nst^2,1);mat2vec_v000(zeros(nst));Cinp(:);mat2vec_v000(Rinp)],'double');
fwrite(FSMO_COV,[mat2vec_v000(Cinp*Psmo*Cinp');mat2vec_v000(Psmo)],'double');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%filtered solution (just for comparison)
Csys=T*[Cti Ctv];
Rsys=(T*Dti)*(T*Dti)';

y=T*sendat(3:end,end);
K=Pfilt*Csys'*inv(Csys*Pfilt*Csys'+Rsys);
Pfilt=(eye(nst)-K*Csys)*Pfilt;
xfilt=xfilt+K*(y-Csys*xfilt);
ufilt=MLS*sendat(3:end,in)-Cinp*xfilt;

fwrite(FFIL_RES,[xfilt;ufilt],'double');
fwrite(FFIL_ERR,[zeros(nst^2,1);mat2vec_v000(zeros(nst));Cinp(:);mat2vec_v000(Rinp)],'double');
fwrite(FFIL_COV,[mat2vec_v000(Cinp*Pfilt*Cinp');mat2vec_v000(Pfilt)],'double');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fclose all;

% return;
%%%%% plot results
sendat=readbin_v000([DIRNAME '\imu.bin'],6);
mlsdat=MLS*sendat(3:6,:);
smodat=readbin_v000([DIRNAME '\smores.bin'],18);
fildat=readbin_v000([DIRNAME '\filres.bin'],18);
smocov=readbin_v000([DIRNAME '\smocov.bin'],nst*(nst+1)/2+3);
filcov=readbin_v000([DIRNAME '\filcov.bin'],nst*(nst+1)/2+3);

figure(1);
plot(mlsdat(1,:),'k')
hold on
plot(fildat(17,:))
plot(smodat(17,:),'r')
plot(filcov(1,:).^0.5,'--')
plot(smocov(1,:).^0.5,'r--')

figure(2);
plot(mlsdat(2,:),'k')
hold on
plot(fildat(18,:))
plot(smodat(18,:),'r')
plot(filcov(2,:).^0.5,'--')
plot(smocov(2,:).^0.5,'r--')