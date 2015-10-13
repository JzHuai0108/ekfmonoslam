function [SenErrDef, IniErrDef, ObsErrDef, ini_pva, wander, dt]=sys_path(dir_name,record)
rstat=[134 26 345 476 5456;654 7345 8345 9678 140];       %control how random the errors are generated
%rstat=rand(2,5)*30000;
dt=1/100;
D2R=pi/180;

%%%imu error parameters (discrete time param. defined at 100Hz)
%Accelerometers (defined in m-sec^2)
SenErrDef(1).A=1;
SenErrDef(1).B=3e-006;
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

% SenErrDef(6).A=1;
% SenErrDef(6).B=0;
% SenErrDef(6).C=1;
% SenErrDef(6).D=0;
% SenErrDef(6).sP=0;
% SenErrDef(1)=SenErrDef(6);
% SenErrDef(2)=SenErrDef(6);
% SenErrDef(3)=SenErrDef(6);
% SenErrDef(4)=SenErrDef(6);
% SenErrDef(5)=SenErrDef(6);

%%initial err defs
IniErrDef.pos_sP=diag([1/6e6,1/6e6,1]);
IniErrDef.vel_sP=eye(3)*0.01;
IniErrDef.att_sP=eye(3)*(10*pi/180);
IniErrDef.wander_sP=eye(2)*2;

%3D simulation in NED
%initial nav values
vel_b=zeros(3,1);
att_n=[-3;5;110]*D2R;
pos_n=[39.88864*D2R;32.78002*D2R;1102]; %51.081755-114.1360749

%%observation err defs
[Rn, Re, g, sL, cL, WIE_E]=geoparam_v000(pos_n);
Sn=Rn+pos_n(3);
Se=(Re+pos_n(3))*cos(pos_n(1));
ObsErrDef.sR=diag([0.5/Sn, 0.5/Se, 1.5]); %position error std (m)


%%%%%%%%%%% Generate path %%%%%%%%%%%%%%%%
if (record)
    %%motion definitions
    mot_def(1,:)=[1 0 0 0 0 5];
    mot_def(2,:)=[5 0 0 0 2 10];
    mot_def(3,:)=[1 0 0 0 0 5];
    mot_def(4,:)=[3 0 0 90*D2R 0 15];
    mot_def(5,:)=[1 0 0 0 0 5];
    mot_def(6,:)=[3 0 0 90*D2R 0 15];
    mot_def(7,:)=[1 0 0 0 0 5];
    mot_def(8,:)=[3 0 0 -90*D2R 0 15];
    mot_def(9,:)=[1 0 0 0 0 10];
    mot_def(10,:)=[3 0 0 -90*D2R 0 15];
    mot_def(11,:)=[1 0 0 0 0 5];
    mot_def(12,:)=[5 0 0 0 0 10];

    PathGen_v002(dir_name, [pos_n, vel_b, att_n], mot_def, [1 1/dt;1 0.1], 0, [], [], [], rstat(:,1));


    % %%%Add error
    % %Imu errors
    AddIMUErr_v000([dir_name 'mimu.bin'], [dir_name 'imu.bin'], [dir_name 'imuerr.bin'], 7, 2:7, SenErrDef, [], rstat(:,2));
    AddObsErr_v000([dir_name 'gps.bin'], 7, [dir_name 'obs_pos.bin'], [5 6 7], ObsErrDef, rstat(:,3));

    %AddIMUErr_v000([dir_name 'mimu.bin'], [dir_name 'imu.bin'], [dir_name 'imuerr.bin'], 7, 2:7, [], [], rstat(:,2));
    %AddObsErr_v000([dir_name 'gps.bin'], 7, [dir_name 'obs_pos.bin'], [5 6 7], [], rstat(:,3));
end

%initialization errors
randn('state',rstat(:,4));
pos_err=IniErrDef.pos_sP*randn(3,1);
vel_err=IniErrDef.vel_sP*randn(3,1);
att_err=IniErrDef.att_sP*randn(3,1);
Cerr=euler2dcm_v000(att_err);
Cbg=euler2dcm_v000(att_n);
Cbg_err=Cerr'*Cbg;
att_ini=dcm2euler_v000(Cbg_err);

ini_pva=[pos_n+pos_err vel_b+vel_err att_ini];  %%Attitude will be initialized using acc data
%wander=[sin(att_ini(3));cos(att_ini(3))];
wander=[sin(105*pi/180);cos(105*pi/180)];
