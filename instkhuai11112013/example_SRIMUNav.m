clear;
PATH='C:\Users\yigiter\Desktop\Dropbox\SRIMU\Fusion\data';
D2R=pi/180;

%%% 2d multi imu simulation
dt=1/25;

%initial nav values
pos_n=[0.891800000000;-1.991979000000;1112.000000000000];
[Rn, Re, g, sL, cL, WIE_E]=geoparam_v000(pos_n);
att_n=[0;0;60]*D2R;
Cbn=euler2dcm_v000(att_n);
vel_b=[0;0;0];
pos_dn=[0;0;0];

%Sensor configurations
sen_conf(1,:)=[0 0 0 cos(0) sin(0) 0 1];
sen_conf(2,:)=[0 0 0 cos(pi/8) sin(pi/8) 0 1];
sen_conf(3,:)=[0 0 0 cos(pi/4) sin(pi/4) 0 1];
sen_conf(4,:)=[0 0 0 cos(-pi/8) sin(-pi/8) 0 1];
sen_conf(5,:)=[0 0 0 cos(-pi/4) sin(-pi/4) 0 1];

sen_conf(6,:)=[0 0 0 cos(pi) sin(pi) 0 1];
sen_conf(7,:)=[0 0 0 cos(pi-pi/8) sin(pi-pi/8) 0 1];
sen_conf(8,:)=[0 0 0 cos(pi-pi/4) sin(pi-pi/4) 0 1];

sen_conf(9,:)=[0 0 0 0 0 -1 2];
sen_conf(10,:)=[0 0 0 0 0 1 2];
sen_conf(11,:)=[0 0 0 0 0 1 2];
act_sen=1:11;
nsen=size(sen_conf,1);


%IMU stability error definitions
LoadErrorDefs2();
for (in=1:length(sen_errdef))   %convert error models according to the new sensor frequency
    mx_a=sen_errdef(in).B*sen_errdef(in).B';
    mx_b=sen_errdef(in).D*sen_errdef(in).D';
    [sen_errdef(in).A, mx_c, mx_d]=dc2dc_v000(sen_errdef(in).A, mx_a, mx_b, 1/100, 2, dt);
    sen_errdef(in).B=mx_c.^0.5;
    sen_errdef(in).D=mx_d.^0.5;
end

% %%%%%%%%%%% Generate path %%%%%%%%%%%%%%%%
%motion definitions
mot_def(1,:)=[1 0 0 0 0 1];
mot_def(2,:)=[5 0 0 0 2 10];
mot_def(3,:)=[1 0 0 0 0 10];
mot_def(4,:)=[3 0 0 90*D2R 0 15];
mot_def(5,:)=[1 0 0 0 0 5];
mot_def(6,:)=[3 0 0 90*D2R 0 15];
mot_def(7,:)=[1 0 0 0 0 5];
mot_def(8,:)=[5 0 0 0 0 10];

%sensor/imu configurations
vr_a=dir2euler_v000(sen_conf(:,4:6));
xconf=[zeros(nsen,3) vr_a, sen_conf(:,7)];

%generate trajectory
load([PATH '\rstat2.mat']);
% rstat=rstat*0;
rstat(:,1)=PathGen_v001(PATH, [pos_n, vel_b, att_n], mot_def, [1 1/dt;1 0.025], 1, xconf, [], [], rstat(:,1));

% %Add error
rstat(:,2)=AddIMUErr_v000([PATH '\mimu.bin'], [PATH '\imu.bin'], [PATH '\imuerr.bin'], nsen+7, [8:nsen+7], sen_errdef(act_sen), tem_mod, rstat(:,2));
rstat(:,4)=AddObsErr_v000([PATH '\gps.bin'], 7, [PATH '\obs.bin'], [5 6 7], obs_errdef(1), rstat(:,4));
[pos_dn, vel_b, Cfbn, rstat(:,4)]=AddIniErr_v000(pos_dn, vel_b, Cbn, ini_errdef, rstat(:,4));
% save([PATH '\rstat2.mat'],'rstat');
% %%%%%%%%%%%%%%%%%%%%%%%%%%


%%read data
sendat=readbin_v000([PATH '\imu.bin'],nsen+2);
obsdat=readbin_v000([PATH '\obs.bin'],4);
trunavdat=readbin_v000([PATH '\mnav.bin'],10);
% trusendat=readbin_v000([PATH '\mimu.bin'],nsen+7);

nst_nav=5;
[nst_sen nst_senTI nst_senTV]=st_size_v000(sen_errdef(act_sen));
trusenerr=readbin_v000([PATH '\imuerr.bin'],nst_sen+1);


%%%%%%%%%%% Start the main part %%%%%%%%%%%%%%%%
%%Subtype models
%%Subtype configurations
M_s1=cp_M_v000(sen_conf(1:5,:));
M_s2=cp_M_v000(sen_conf(6:8,:));
M_s3=cp_M_v000(sen_conf(9:11,:));

%subtype models
[AimuTI_s1, BimuTI_s1, CimuTI_s1, DimuTI_s1, sP0imuTI_s1]=imu_modTI_v000(sen_errdef(act_sen(1)));
[AimuTI_s2, BimuTI_s2, CimuTI_s2, DimuTI_s2, sP0imuTI_s2]=imu_modTI_v000(sen_errdef(act_sen(7)));
[AimuTI_s3, BimuTI_s3, CimuTI_s3, DimuTI_s3, sP0imuTI_s3]=imu_modTI_v000(sen_errdef(act_sen(10)));

[Mls_s1, AimuTI_s1, BimuTI_s1, CimuTI_s1, DimuTI_s1, sP0imuTI_s1]=cp_redsys_v000(M_s1, AimuTI_s1, BimuTI_s1, CimuTI_s1, DimuTI_s1, sP0imuTI_s1);
[Mls_s2, AimuTI_s2, BimuTI_s2, CimuTI_s2, DimuTI_s2, sP0imuTI_s2]=cp_redsys_v000(M_s2, AimuTI_s2, BimuTI_s2, CimuTI_s2, DimuTI_s2, sP0imuTI_s2);
[Mls_s3, AimuTI_s3, BimuTI_s3, CimuTI_s3, DimuTI_s3, sP0imuTI_s3]=cp_redsys_v000(M_s3, AimuTI_s3, BimuTI_s3, CimuTI_s3, DimuTI_s3, sP0imuTI_s3);

AimuTIEq=diagmat_v000(AimuTI_s2,AimuTI_s1,[]);AimuTIEq=diagmat_v000(AimuTI_s3,AimuTIEq,[]);
BimuTIEq=diagmat_v000(BimuTI_s2,BimuTI_s1,[]);BimuTIEq=diagmat_v000(BimuTI_s3,BimuTIEq,[]);
CimuTIEq=diagmat_v000(CimuTI_s2,CimuTI_s1,[]);CimuTIEq=diagmat_v000(CimuTI_s3,CimuTIEq,[]);
DimuTIEq=diagmat_v000(DimuTI_s2,DimuTI_s1,[]);DimuTIEq=diagmat_v000(DimuTI_s3,DimuTIEq,[]);
sP0imuTIEq=diagmat_v000(sP0imuTI_s2,sP0imuTI_s1,[]);sP0imuTIEq=diagmat_v000(sP0imuTI_s3,sP0imuTIEq,[]);

nst_senEq=size(AimuTIEq,1);
MEq=[1 0 0;0 1 0;1 0 0;0 1 0;0 0 1];
[T1Eq MlsEq]=cp_Tparam_v000(MEq,DimuTIEq*DimuTIEq');

Cimu_obsEq=T1Eq*CimuTIEq;
Himu_obsEq=[zeros(size(T1Eq,1),nst_nav) Cimu_obsEq];
Cimu_inpEq=MlsEq*CimuTIEq;
Rimu_obsEq=T1Eq*(DimuTIEq*DimuTIEq')*T1Eq';
Rimu_inpEq=MlsEq*(DimuTIEq*DimuTIEq')*MlsEq';


QimuEq=BimuTIEq*BimuTIEq';
PimuEq=sP0imuTIEq*sP0imuTIEq';
AimuEq=AimuTIEq;


%%%%%%%%%%%%%%%%% Overall mnodel
[AimuTI, BimuTI, CimuTI, DimuTI, sP0imuTI]=imu_modTI_v000(sen_errdef(act_sen));
%sensor obs matrix

M=cp_M_v000(sen_conf(act_sen,:));
[T1 Mls]=cp_Tparam_v000(M,DimuTI*DimuTI');
% [mx_a mx_b]=cp_Tparam_v000(M(1:9,1:2),DimuTI(1:9,1:9)*DimuTI(1:9,1:9)');
% [mx_c mx_d]=cp_Tparam_v000(M(10:11,3),DimuTI(10:11,10:11)*DimuTI(10:11,10:11)');
% 
% mx_e=diagmat_v000(mx_c,mx_a,1);
% T1=mx_e;

%Transformed imu observation matrices
Cimu_obs=T1*CimuTI;
Himu_obs=[zeros(size(T1,1),nst_nav) Cimu_obs];
Cimu_inp=Mls*CimuTI;
Rimu_obs=T1*(DimuTI*DimuTI')*T1';
Rimu_inp=Mls*(DimuTI*DimuTI')*Mls';

Qimu=BimuTI*BimuTI';
Pimu=sP0imuTI*sP0imuTI';
Aimu=AimuTI;


%initial values
pos_ml=pos_dn;
pos_op=pos_dn;
pos_po=pos_dn;
pos_opeq=pos_dn;

vel_ml=vel_b;
vel_op=vel_b;
vel_po=vel_b;
vel_opeq=vel_b;

Cbn_ml=Cbn;
Cbn_op=Cbn;
Cbn_po=Cbn;
Cbn_opeq=Cbn;

%imu err states
ximu_ml=zeros(nst_senEq,1);
ximu_op=zeros(nst_sen,1);
ximu_opeq=zeros(nst_senEq,1);
ximu_poa=zeros(nst_sen,1);
ximu_pob=zeros(nst_sen,1);


%covariance values
Pnav=zeros(nst_nav);
Pnav(1:2,1:2)=ini_errdef.pos_sP(1:2,1:2)*ini_errdef.pos_sP(1:2,1:2)';
Pnav(3:4,3:4)=ini_errdef.vel_sP(1:2,1:2)*ini_errdef.vel_sP(1:2,1:2)';
Pnav(5,5)=ini_errdef.att_sP(3,3)*ini_errdef.att_sP(3,3)';

Pml=zeros(nst_nav+nst_senEq);
Pml(1:nst_nav,1:nst_nav)=Pnav;
Pml(nst_nav+[1:nst_senEq],nst_nav+[1:nst_senEq])=PimuEq;
Pop=zeros(nst_nav+nst_sen);
Pop(1:nst_nav,1:nst_nav)=Pnav;
Pop(nst_nav+[1:nst_sen],nst_nav+[1:nst_sen])=Pimu;
Popeq=zeros(nst_nav+nst_senEq);
Popeq(1:nst_nav,1:nst_nav)=Pnav;
Popeq(nst_nav+[1:nst_senEq],nst_nav+[1:nst_senEq])=PimuEq;
Ppo=zeros(nst_nav+nst_sen);
Ppo(1:nst_nav,1:nst_nav)=Pnav;
Ppo(nst_nav+[1:nst_sen],nst_nav+[1:nst_sen])=Pimu;

%%external observation models
obs_ctr=1;
RExt_obs=obs_errdef(1).sR*obs_errdef(1).sR';
RExt_obs=RExt_obs(1:2,1:2);
HExt_obs=[[1 0 0 0 0;0 1 0 0 0] zeros(2,nst_sen)];
HExt_obsEq=[[1 0 0 0 0;0 1 0 0 0] zeros(2,nst_senEq)];

%%%Start nav
%The first imu outputs and imu observation
imu_eq_s1=Mls_s1*sendat(3:7,2);
imu_eq_s2=Mls_s2*sendat(8:10,2);
imu_eq_s3=Mls_s3*sendat(11:13,2);
sendatEq=[imu_eq_s1;imu_eq_s2;imu_eq_s3];
sendatAll=sendat(3:end,2);
imu_ml=MlsEq*sendatEq;

imu_op=Mls*sendatAll;
imu_opeq=MlsEq*sendatEq;
imu_po=Mls*sendatAll;

imu_obs_po=(T1*sendatAll)-(Cimu_obs*ximu_poa);
imu_obs_op=(T1*sendatAll)-(Cimu_obs*ximu_op);
imu_obs_eq=(T1Eq*sendatEq)-(Cimu_obsEq*ximu_opeq);


%imu kalman
Kimu=Pimu*Cimu_obs'*inv(Cimu_obs*Pimu*Cimu_obs'+Rimu_obs);
Pimu=Pimu-Kimu*Cimu_obs*Pimu;
dz=Kimu*(imu_obs_po);
ximu_poa=ximu_poa+dz;
imu_po=imu_po-Cimu_inp*ximu_poa;

%po ini cov
Ppo(nst_nav+[1:nst_sen],nst_nav+[1:nst_sen])=Pimu;

%op kalman
Kop=Pop*Himu_obs'*inv(Himu_obs*Pop*Himu_obs'+Rimu_obs);
Pop=Pop-Kop*Himu_obs*Pop;
dz=Kop*(imu_obs_op);
pos_op(1:2)=pos_op(1:2)-dz(1:2);
vel_op(1:2)=vel_op(1:2)-dz(3:4);
Cbn_op=(eye(3)+skew([0;0;dz(5)]))*Cbn_op;
ximu_op=ximu_op+dz(nst_nav+[1:nst_sen]);
imu_op=imu_op-Cimu_inp*ximu_op;


%opeq kalman
Kopeq=Popeq*Himu_obsEq'*inv(Himu_obsEq*Popeq*Himu_obsEq'+Rimu_obsEq);
Popeq=Popeq-Kopeq*Himu_obsEq*Popeq;
dz=Kopeq*(imu_obs_eq);
pos_opeq(1:2)=pos_opeq(1:2)-dz(1:2);
vel_opeq(1:2)=vel_opeq(1:2)-dz(3:4);
Cbn_opeq=(eye(3)+skew([0;0;dz(5)]))*Cbn_opeq;
ximu_opeq=ximu_opeq+dz(nst_nav+[1:nst_senEq]);
imu_opeq=imu_opeq-Cimu_inpEq*ximu_opeq;

%%%%
%%Apply additional obs here
%%%%

%%debug
%%%%debug variables
debctr=0;
debimu=zeros(12,size(sendat,2));
debnav=zeros(20,size(sendat,2));
deberr_ml=zeros(nst_senEq+nst_nav,size(sendat,2));
deberr_op=zeros(nst_sen+nst_nav,size(sendat,2));
deberr_po=zeros(nst_sen+nst_nav,size(sendat,2));
deberr_eq=zeros(nst_senEq+nst_nav,size(sendat,2));
debP_ml=zeros(nst_senEq+nst_nav,size(sendat,2));
debP_op=zeros(nst_sen+nst_nav,size(sendat,2));
debP_opeq=zeros(nst_senEq+nst_nav,size(sendat,2));
debP_po=zeros(nst_sen+nst_nav,size(sendat,2));

debctr=debctr+1;
yy=1;
debimu(:,debctr)=[imu_ml;imu_op;imu_po;imu_opeq];
debnav(:,debctr)=[pos_ml(1:2);vel_ml(1:2);dcm2heading_v000(Cbn_ml);pos_op(1:2);vel_op(1:2);dcm2heading_v000(Cbn_op);pos_po(1:2);vel_po(1:2);dcm2heading_v000(Cbn_po);pos_opeq(1:2);vel_opeq(1:2);dcm2heading_v000(Cbn_opeq)];
vr_a=pos_ml(1:2)-trunavdat(2:3,yy);
vr_b=vel_ml(1:2)-trunavdat(5:6,yy);
vr_c=dcm2heading_v000(Cbn_ml'*euler2dcm_v000(trunavdat(8:10,yy)));
vr_d=ximu_ml;%-trusenerr(2:end,yy);
deberr_ml(:,debctr)=[vr_a;vr_b;vr_c;vr_d];
debP_ml(:,debctr)=diag(Pml).^0.5;
vr_a=pos_op(1:2)-trunavdat(2:3,yy);
vr_b=vel_op(1:2)-trunavdat(5:6,yy);
vr_c=dcm2heading_v000(Cbn_op'*euler2dcm_v000(trunavdat(8:10,yy)));
vr_d=ximu_op;%-trusenerr(2:end,yy);
deberr_op(:,debctr)=[vr_a;vr_b;vr_c;vr_d];
debP_op(:,debctr)=diag(Pop).^0.5;
vr_a=pos_po(1:2)-trunavdat(2:3,yy);
vr_b=vel_po(1:2)-trunavdat(5:6,yy);
vr_c=dcm2heading_v000(Cbn_po'*euler2dcm_v000(trunavdat(8:10,yy)));
vr_d=ximu_pob-trusenerr(2:end,yy);
deberr_po(:,debctr)=[vr_a;vr_b;vr_c;vr_d];
debP_po(:,debctr)=diag(Ppo).^0.5;
vr_a=pos_opeq(1:2)-trunavdat(2:3,yy);
vr_b=vel_opeq(1:2)-trunavdat(5:6,yy);
vr_c=dcm2heading_v000(Cbn_opeq'*euler2dcm_v000(trunavdat(8:10,yy)));
vr_d=ximu_opeq;%-trusenerr(2:end,yy);
deberr_opeq(:,debctr)=[vr_a;vr_b;vr_c;vr_d];
debP_opeq(:,debctr)=diag(Popeq).^0.5;

%%%%%%%%%%%%%%

for in=3:size(sendat,2)
    %run the strapdowns
    [Cbn_ml, vel_ml, pos_ml]=strapdown_pln_dcm_v000(Cbn_ml, vel_ml, pos_ml, [imu_ml(1:2);-g], [0;0;imu_ml(3)], g, dt, 1);
    [Cbn_op, vel_op, pos_op]=strapdown_pln_dcm_v000(Cbn_op, vel_op, pos_op, [imu_op(1:2);-g], [0;0;imu_op(3)], g, dt, 1);
    [Cbn_opeq, vel_opeq, pos_opeq]=strapdown_pln_dcm_v000(Cbn_opeq, vel_opeq, pos_opeq, [imu_opeq(1:2);-g], [0;0;imu_opeq(3)], g, dt, 1);
    [Cbn_po, vel_po, pos_po]=strapdown_pln_dcm_v000(Cbn_po, vel_po, pos_po, [imu_po(1:2);-g], [0;0;imu_po(3)], g, dt, 1);    
    
    %%%new imu interrupt at this point
    sendatAll=sendat(3:end,in);
    imu_eq_s1=Mls_s1*sendat(3:7,in);
    imu_eq_s2=Mls_s2*sendat(8:10,in);
    imu_eq_s3=Mls_s3*sendat(11:13,in);
    sendatEq=[imu_eq_s1;imu_eq_s2;imu_eq_s3];
    
    
    %%%imu KF  
    %Predict
    Pimu=Aimu*Pimu*Aimu+Qimu';
    ximu_poa=Aimu*ximu_poa;
    %update
    imu_obs_po=(T1*sendatAll)-(Cimu_obs*ximu_poa);
    Kimu=Pimu*Cimu_obs'*inv(Cimu_obs*Pimu*Cimu_obs'+Rimu_obs);
    Pimu=Pimu-Kimu*Cimu_obs*Pimu;
    dz=Kimu*(imu_obs_po);
    ximu_poa=ximu_poa+dz;
        
    %form the Nav model
    [Anav N]=sys_pln_dcm_v000(pos_ml, vel_ml, Cbn_ml, [imu_ml(1:2);-g], [0 0 imu_ml(3)], g, 1, 1, dt);
    Anav([3 6 7 8],:)=[];
    Anav(:,[3 6 7 8])=[];
    N([3 6 7 8],:)=[];
    N(:, [3 4 5])=[];
      
    %Form the new imu model for po
    Aimu_po=Aimu-Kimu*Cimu_obs*Aimu;
    Qimu_po=Kimu*Rimu_obs*Kimu'+ (eye(nst_sen)-Kimu*Cimu_obs)*Qimu*(eye(nst_sen)-Kimu*Cimu_obs)';
    
    %%%%KF predictions
    %overal system matrices
    Aml=[Anav N*Cimu_inpEq*dt; zeros(nst_senEq, nst_nav) AimuEq];
    Qml=diagmat_v000(QimuEq, N*Rimu_inpEq*dt*N'*dt,1);

    Aop=[Anav N*Cimu_inp*dt; zeros(nst_sen, nst_nav) Aimu];
    Qop=diagmat_v000(Qimu, N*Rimu_inp*dt*N'*dt,1);
    
    Aopeq=Aml;
    Qopeq=Qml;
        
    Apo=[Anav N*Cimu_inp*dt; zeros(nst_sen, nst_nav) Aimu_po];
    Qpo=diagmat_v000(Qimu_po, N*Rimu_inp*dt*N'*dt,1);
    
    %Predictions
    ximu_ml=AimuEq*ximu_ml;
    ximu_op=Aimu*ximu_op;
    ximu_opeq=AimuEq*ximu_opeq;
    ximu_pob=Aimu_po*ximu_pob;
    Pml=Aml*Pml*Aml'+Qml;
    Pop=Aop*Pop*Aop'+Qop;
    Popeq=Aopeq*Popeq*Aopeq'+Qopeq;
    Ppo=Apo*Ppo*Apo'+Qpo;
    
    %op kalman  
    imu_obs_op=(T1*sendatAll)-(Cimu_obs*ximu_op);
    Kop=Pop*Himu_obs'*inv(Himu_obs*Pop*Himu_obs'+Rimu_obs);
    Pop=Pop-Kop*Himu_obs*Pop;
    dz=Kop*(imu_obs_op);
    pos_op(1:2)=pos_op(1:2)-dz(1:2);
    vel_op(1:2)=vel_op(1:2)-dz(3:4);
    Cbn_op=(eye(3)+skew([0;0;dz(5)]))*Cbn_op;
    ximu_op=ximu_op+dz(nst_nav+[1:nst_sen]);
    
    
    %opeq kalman  
    imu_obs_eq=(T1Eq*sendatEq)-(Cimu_obsEq*ximu_opeq);
    Kopeq=Popeq*Himu_obsEq'*inv(Himu_obsEq*Popeq*Himu_obsEq'+Rimu_obsEq);
    Popeq=Popeq-Kopeq*Himu_obsEq*Popeq;
    dz=Kopeq*(imu_obs_eq);
    pos_opeq(1:2)=pos_opeq(1:2)-dz(1:2);
    vel_opeq(1:2)=vel_opeq(1:2)-dz(3:4);
    Cbn_opeq=(eye(3)+skew([0;0;dz(5)]))*Cbn_opeq;
    ximu_opeq=ximu_opeq+dz(nst_nav+[1:nst_senEq]);

    %new imu outputs
    imu_op=Mls*sendatAll-Cimu_inp*(ximu_op);
    imu_opeq=MlsEq*sendatEq-Cimu_inpEq*(ximu_opeq);
    imu_ml=MlsEq*sendatEq-Cimu_inpEq*(ximu_ml);
    imu_po=Mls*sendatAll-Cimu_inp*(ximu_poa+ximu_pob);
    
    
    %%%%
    %%Apply additional obs here
    if (obs_ctr<=size(obsdat,2))
        if (obsdat(1,obs_ctr)==(in-2)) %-2 is just because of my indexing convention
            pos_obs=obsdat(2:4,obs_ctr);

            %ml kalman
            pos_obs_ml=pos_ml-pos_obs;
            Kml=Pml*HExt_obsEq'*inv(HExt_obsEq*Pml*HExt_obsEq'+RExt_obs);
            Pml=Pml-Kml*HExt_obsEq*Pml;
            dz=Kml*(pos_obs_ml(1:2));
            pos_ml(1:2)=pos_ml(1:2)-dz(1:2);
            vel_ml(1:2)=vel_ml(1:2)-dz(3:4);
            Cbn_ml=(eye(3)+skew([0;0;dz(5)]))*Cbn_ml;
            ximu_ml=ximu_ml+dz(nst_nav+[1:nst_senEq]);

            %op kalman  
            pos_obs_op=pos_op-pos_obs;
            Kop=Pop*HExt_obs'*inv(HExt_obs*Pop*HExt_obs'+RExt_obs);
            Pop=Pop-Kop*HExt_obs*Pop;
            dz=Kop*(pos_obs_op(1:2));
            pos_op(1:2)=pos_op(1:2)-dz(1:2);
            vel_op(1:2)=vel_op(1:2)-dz(3:4);
            Cbn_op=(eye(3)+skew([0;0;dz(5)]))*Cbn_op;
            ximu_op=ximu_op+dz(nst_nav+[1:nst_sen]);
            
            
            %opeq kalman  
            pos_obs_eq=pos_opeq-pos_obs;
            Kopeq=Popeq*HExt_obsEq'*inv(HExt_obsEq*Popeq*HExt_obsEq'+RExt_obs);
            Popeq=Popeq-Kopeq*HExt_obsEq*Popeq;
            dz=Kopeq*(pos_obs_eq(1:2));
            pos_opeq(1:2)=pos_opeq(1:2)-dz(1:2);
            vel_opeq(1:2)=vel_opeq(1:2)-dz(3:4);
            Cbn_opeq=(eye(3)+skew([0;0;dz(5)]))*Cbn_opeq;
            ximu_opeq=ximu_opeq+dz(nst_nav+[1:nst_senEq]);

            %po kalman  
            pos_obs_po=pos_po-pos_obs;
            Kpo=Ppo*HExt_obs'*inv(HExt_obs*Ppo*HExt_obs'+RExt_obs);
            Ppo=Ppo-Kpo*HExt_obs*Ppo;
            dz=Kpo*(pos_obs_po(1:2));
            pos_po(1:2)=pos_po(1:2)-dz(1:2);
            vel_po(1:2)=vel_po(1:2)-dz(3:4);
            Cbn_po=(eye(3)+skew([0;0;dz(5)]))*Cbn_po;
            ximu_pob=ximu_pob+dz(nst_nav+[1:nst_sen]);

            obs_ctr=obs_ctr+1;
        end
    end
    %%%%
    
    %%debug
    debctr=debctr+1;
    yy=in;
    debimu(:,debctr)=[imu_ml;imu_op;imu_po;imu_opeq];
    debnav(:,debctr)=[pos_ml(1:2);vel_ml(1:2);dcm2heading_v000(Cbn_ml);pos_op(1:2);vel_op(1:2);dcm2heading_v000(Cbn_op);pos_po(1:2);vel_po(1:2);dcm2heading_v000(Cbn_po);pos_opeq(1:2);vel_opeq(1:2);dcm2heading_v000(Cbn_opeq)];
    vr_a=pos_ml(1:2)-trunavdat(2:3,yy-1);
    vr_b=vel_ml(1:2)-trunavdat(5:6,yy-1);
    vr_c=dcm2heading_v000(Cbn_ml'*euler2dcm_v000(trunavdat(8:10,yy-1)));
    vr_d=ximu_ml;%-trusenerr(2:end,yy);
    deberr_ml(:,debctr)=[vr_a;vr_b;vr_c;vr_d];
    debP_ml(:,debctr)=diag(Pml).^0.5;
    vr_a=pos_op(1:2)-trunavdat(2:3,yy-1);
    vr_b=vel_op(1:2)-trunavdat(5:6,yy-1);
    vr_c=dcm2heading_v000(Cbn_op'*euler2dcm_v000(trunavdat(8:10,yy-1)));
    vr_d=ximu_op;%-trusenerr(2:end,yy);
    deberr_op(:,debctr)=[vr_a;vr_b;vr_c;vr_d];
    debP_op(:,debctr)=diag(Pop).^0.5;
    vr_a=pos_po(1:2)-trunavdat(2:3,yy-1);
    vr_b=vel_po(1:2)-trunavdat(5:6,yy-1);
    vr_c=dcm2heading_v000(Cbn_po'*euler2dcm_v000(trunavdat(8:10,yy-1)));
    vr_d=ximu_pob-trusenerr(2:end,yy);
    deberr_po(:,debctr)=[vr_a;vr_b;vr_c;vr_d];
    debP_po(:,debctr)=diag(Ppo).^0.5;
    vr_a=pos_opeq(1:2)-trunavdat(2:3,yy-1);
    vr_b=vel_opeq(1:2)-trunavdat(5:6,yy-1);
    vr_c=dcm2heading_v000(Cbn_opeq'*euler2dcm_v000(trunavdat(8:10,yy-1)));
    vr_d=ximu_opeq;%-trusenerr(2:end,yy);
    deberr_opeq(:,debctr)=[vr_a;vr_b;vr_c;vr_d];
    debP_opeq(:,debctr)=diag(Popeq).^0.5;
end
fclose all;

% pos_ml
% pos_op
% pos_po

% vel_ml
% vel_op
% vel_po

% dcm2euler_v000(Cbn_ml)
% dcm2euler_v000(Cbn_op)
% dcm2euler_v000(Cbn_po)


figure(13);
plot(trunavdat(2,:),trunavdat(3,:));
hold on
plot(debnav(1,1:debctr),debnav(2,1:debctr),'r');
plot(debnav(6,1:debctr),debnav(7,1:debctr),'k');
plot(debnav(11,1:debctr),debnav(12,1:debctr),'g');
plot(debnav(16,1:debctr),debnav(17,1:debctr),'k');
grid;
xlabel('X Pos. (m)');
ylabel('Y Pos. (m)');
legend('True Trajectory', 'x^{wls}', 'x^{op}', 'x^{po}', 'x^{rop}')


figure (14);
tt=[1:debctr]/25;
plot(tt,deberr_op(4,1:debctr));
hold on
plot(tt,deberr_ml(4,1:debctr),'r');
plot(tt,deberr_po(4,1:debctr),'k');
plot(tt,deberr_opeq(4,1:debctr));
plot(tt,debP_op(4,1:debctr),'--');
plot(tt,debP_ml(4,1:debctr),'--r');
plot(tt,debP_po(4,1:debctr),'--k');
plot(tt,debP_opeq(4,1:debctr),'--');
grid;
xlabel('Time (Sec)');
ylabel('Y Vel. Error (m/sec)');
legend('True Trajectory', 'x^{wls}', 'x^{op}', 'x^{po}', 'x^{rop}')
legend('x^{op}', 'x^{wls}', 'x^{po}', 'x^{rop}','StdDev of x^{op}','StdDev of x^{wls}','StdDev of x^{po}')


return;

figure;
tt=[1:debctr]/25;
plot(tt,trunavdat(5,1:debctr));
hold on
plot(tt,debnav(3,1:debctr),'r');
plot(tt,debnav(8,1:debctr),'k');
plot(tt,debnav(13,1:debctr),'g');
plot(tt,debnav(18,1:debctr),'k');
grid;


figure;

sr_a=1;
plot(tt, deberr_op(5+sr_a,1:debctr))
hold on;
plot([1:size(trusenerr,2)]/25, trusenerr(sr_a+1,:),'r');
plot(tt,debP_po(5+sr_a,1:debctr),'--');
grid;
xlabel('Time (sec)');
ylabel('Acc_1 Bias (m/sec^2)');
legend('OP Sol', 'True Bias', 'KF Std. Dev.');

vr_a=Cimu_inp*deberr_op(6:end,1:debctr);
vr_b=Cimu_inpEq*deberr_ml(6:end,1:debctr);
vr_c=Cimu_inp*trusenerr(2:end,1:debctr);
figure;
plot(tt,vr_a(3,:));
hold on;
plot(tt,vr_b(3,:));
plot(tt,vr_c(3,:),'r');
grid;
legend('Op Sol.','Reduced Order ML.','True Eff Gyro Bias');
xlabel('Time (sec)');
ylabel('Eff. Gyro Bias (rad/sec)');

return;