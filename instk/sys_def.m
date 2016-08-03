function [SenErrDef, M, pva, imu_out, dt]=sys_def(DIRNAME)
dt=1/100;
DT=60; %sec.

%%Randomize
rstr=RandStream.getDefaultStream();
reset(rstr,83474);
rstat=rstr.State;

%Sensor Error Definitions (at 100Hz)
%Accelerometers (defined in m-sec^2)
SenErrDef(1).A=0.999988000072;
SenErrDef(1).B=3.49997900010378e-006;
SenErrDef(1).C=1;
SenErrDef(1).D=1e-3/sqrt(0.01);
SenErrDef(1).sP=sqrt(5.1e-007)*1000;
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

%Extra Sensors
% SenErrDef(7)=SenErrDef(1);
% SenErrDef(8)=SenErrDef(4);

%%Sensor Configuration matrix
% M=[eye(6);0.5 0 0.5 0 0 0;0 0 0 0 0 1];
% M(7,:)=M(7,:)/norm(M(7,:));
M=eye(6);

%pva
pos_n=[0.891800000000;-1.991979000000;1112.000000000000];
att_n=[0;0;90]*pi/180;
vel_b=[0;0;0];
pva=[pos_n,vel_b,att_n];

%Real IMU outputs
[Rn, Re, g, sL, cL, WIE_E]=geoparam_v000(pos_n);
Cbn=euler2dcm_v000(att_n);
acc=Cbn'*[0;0;-g];
gyro=Cbn'*[WIE_E*cL;0;-WIE_E*sL];
imu_out=M*[acc;gyro];

%%%Generate sensor outputs
AddIMUErr_v001([], [DIRNAME '\imu.bin'], [DIRNAME '\imuerr.bin'], DT/dt, imu_out, SenErrDef, [], rstat);
