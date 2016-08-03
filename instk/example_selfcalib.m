clear;
DIRNAME='C:\Benisil';

%Define the System and Generate Stationary Sensor Outputs
[SenErrDef, M, pva, imu_th, dt]=sys_def(DIRNAME);

%load Imu outputs
imu_data=readbin_v000([DIRNAME '\imu.bin'],size(M,1)+2);

%Imu model
[Aimu Bimu Cimu Dimu sPimu]=imu_modTI_v000(SenErrDef);

%compute u and self calibrate
%[u, ximu, Pimu, Pu, Pux]=selfcalib1_v000(imu_data(3:end,:), Aimu, Bimu*Bimu',Cimu, Dimu*Dimu', sPimu*sPimu', M, []);
sensor_ref=zeros(size(M,2),4);
sensor_ref(4:6,:)=ones(3,1)*[1 0 0.001 0.0001];
sensor_ref=[];
[u, ximu, Pimu, Pu, Pux]=selfcalib1_v001(imu_data(3:end,:), Aimu, Bimu*Bimu',Cimu, Dimu*Dimu', inv(sPimu*sPimu'), M, sensor_ref);

%[u1, ximu1, Pimu1, Pu1, Pux1]=selfcalib2_v000(imu_data(3:end,:), Aimu, Bimu*Bimu',Cimu, Dimu*Dimu', sPimu*sPimu', M, sensor_ref);
