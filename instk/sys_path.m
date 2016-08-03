function [SenErrDef, IniErrDef, ObsErrDef, ini_pva, dt]=sys_path(dir_name)
rstat=randn(2,5);       %control how random the errors are generated
dt=1/100; %simulation at 100hz
D2R=pi/180;

%%%imu error parameters (discrete time param. defined at 100Hz)
%Accelerometers (defined in m-sec^2)
SenErrDef(1).A=0.999988000072;
SenErrDef(1).B=3.49997900010378e-006;
SenErrDef(1).C=1;
SenErrDef(1).D=1e-3/sqrt(0.01);
SenErrDef(1).sP=sqrt(5.1e-007);
SenErrDef(1).tparam=[];

SenErrDef(2)=SenErrDef(1);
SenErrDef(3)=SenErrDef(1);

%Gyroscopes (defined in rad/sec)
SenErrDef(4).A=0.999980000199999;
SenErrDef(4).B=1.9999800001658e-006;
SenErrDef(4).C=1;
SenErrDef(4).D=0.0006/sqrt(0.01);
SenErrDef(4).sP=sqrt(1e-007);
SenErrDef(4).tparam=[];  %first param=sP of temperature scale factor (RC), second=B of temperature RW

SenErrDef(5)=SenErrDef(4);
SenErrDef(6)=SenErrDef(4);

%%observation err defs
ObsErrDef(1).sR=diag([0.5/6e6,0.5/6e6,1.5]); %position error std (m)
%ObsErrDef(2).sR=eye(3)*0.25; %(zero) velocity error (m/sec)

%%initial err defs
IniErrDef.pos_sP=diag([0.5/6e6,0.5/6e6,1]);
IniErrDef.vel_sP=eye(3)*0.01;
IniErrDef.att_sP=eye(3)*0.01;

%3D simulation in NED
%initial nav values
pos_n=[0.891800000000;-1.991979000000;1112.000000000000];
att_n=[0;0;60]*D2R;
vel_b=[0;0;0];

%%%%%%%%%%% Generate path %%%%%%%%%%%%%%%%
%%motion definitions
mot_def(1,:)=[1 0 0 0 0 1];
mot_def(2,:)=[5 0 0 0 2 10];
mot_def(3,:)=[1 0 0 0 0 10];
mot_def(4,:)=[3 0 0 90*D2R 0 15];
mot_def(5,:)=[1 0 0 0 0 5];
mot_def(6,:)=[3 0 0 90*D2R 0 15];
mot_def(7,:)=[1 0 0 0 0 5];
mot_def(8,:)=[5 0 0 0 0 10];


% %generate trajectory
PathGen_v001(dir_name, [pos_n, vel_b, att_n], mot_def, [1 1/dt;1 1/15], 0, [], [], [], rstat(:,1));

% %%%Add error
% %Imu errors
AddIMUErr_v000([dir_name '\mimu.bin'], [dir_name '\imu.bin'], [dir_name '\imuerr.bin'], 7, 2:7, SenErrDef, [], rstat(:,2));
% %obs errors
AddObsErr_v000([dir_name '\gps.bin'], 7, [dir_name '\obs_pos.bin'], [5 6 7], ObsErrDef, rstat(:,3));

%initialization errors
randn('state',rstat(:,4));
pos_err=IniErrDef.pos_sP*randn(3,1);
vel_err=IniErrDef.vel_sP*randn(3,1);
att_err=IniErrDef.att_sP*randn(3,1);
Cerr=euler2dcm_v000(att_err);
Cbn=euler2dcm_v000(att_n);
Cbn_err=Cerr'*Cbn;
att_ini=dcm2euler_v000(Cbn_err);
ini_pva=[pos_n+pos_err  vel_b+vel_err att_ini];
%ini_pva=[pos_n  vel_b att_n];
