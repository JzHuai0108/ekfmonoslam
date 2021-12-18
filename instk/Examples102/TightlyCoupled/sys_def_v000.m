%%% System Definition Parameters
function [SenErrDef, ClkErrDef, IniErrDef, ObsErrDef, ini_pva, dt]=sys_def()
dt=1/400; %400 Hz

%Initial PVA
Llh=[40*pi/180;35*pi/180;1099];
vel_n=[-0.02;-0.02;0.000];
att_n=[0.0;0.002;1.5];
ini_pva=[Llh vel_n att_n];

%Initial PVA covariance
IniErrDef.pos_sP=eye(3)*10;
IniErrDef.vel_sP=eye(3)*0.2;
IniErrDef.att_sP=eye(3)*0.03;

%%IMU error definitions (in continious time)
%Accelerometers
SenErrDef(1).A=-0.0005;
SenErrDef(1).B=1e-3;
SenErrDef(1).C=1;
SenErrDef(1).D=5e-2;
SenErrDef(1).sP=0.1;
SenErrDef(1).tparam=[];

SenErrDef(2)=SenErrDef(1);
SenErrDef(3)=SenErrDef(1);

%Gyroscopes
SenErrDef(4).A=-0.0003;
SenErrDef(4).B=1e-5;
SenErrDef(4).C=1;
SenErrDef(4).D=1e-2;
SenErrDef(4).sP=1e-2;
SenErrDef(4).tparam=[]; 
SenErrDef(5)=SenErrDef(4);
SenErrDef(6)=SenErrDef(4);

%%Clock Model
ClkErrDef(1).A=[0 1;0 0]; %First state is the bias second state is the drift
ClkErrDef(1).B=diag([0.25,0.5]);
ClkErrDef(1).C=[1 0];
ClkErrDef(1).D=0;
ClkErrDef(1).sP=diag([50000.0, 10]);

%Standard Deviation of PseudoRanges
ObsErrDef.sR=5;
