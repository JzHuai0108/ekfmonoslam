function testAddIMUError(trueImuFile, noisyImuFile, imuType, sampleInterval)
% add VRW, ARW and accelerometer bias noise and gyro bias noise to imu data
% trueImuFile contains ground truth imu samples, i.e. the output of SE3interpolation
% noisyImuFile contains output imu samples with noise
% imuType please refer to instk/imu_err_defs_v000
% sampleInterval, imu data sample interval

addpath('..\..\instk\');
samples=load(trueImuFile);
true_inp=samples(:, 2:7)'; % true input
ErrDefs=cell(6,1);
reducefactor=2;
IMU_ERRDEF=imu_err_defs_v000(imuType);
for i=1:3
  ErrDefs{i}.Cti=zeros(1,12);  
  ErrDefs{i}.Cti(i)=1; 
  ErrDefs{i}.Ctv=zeros(1,12);     
  ErrDefs{i}.R=IMU_ERRDEF.acc_vrw/reducefactor^2; % where D*D'=R, 3dm gx3-35 accel noise density
  ErrDefs{i}.M=zeros(1, 6);
  ErrDefs{i}.M(i)=1;
end
for i=4:6
  ErrDefs{i}.Cti=zeros(1,12); 
  ErrDefs{i}.Cti(i)=1;
  ErrDefs{i}.Ctv=zeros(1,12);     
  ErrDefs{i}.R=IMU_ERRDEF.gyro_arw/reducefactor^2; % where D*D'=R, 3dm gx3-35 gyro noise density
  ErrDefs{i}.M=zeros(1, 6);
  ErrDefs{i}.M(i)=1;
end
StateModel.A=zeros(12);
StateModel.A(1:3,1:3)=-eye(3)/IMU_ERRDEF.acc_bias_Tc;
StateModel.Q=zeros(12);
StateModel.Q(1:3,1:3)=eye(3)*IMU_ERRDEF.acc_bias_Q/reducefactor^2;
StateModel.Q(4:6,4:6)=eye(3)*IMU_ERRDEF.gyro_bias_Q/reducefactor^2;
STD0=zeros(12);
STD0(1:3,1:3)=eye(3)*sqrt(IMU_ERRDEF.acc_bias_var)/reducefactor;
STD0(4:6,4:6)=eye(3)*sqrt(IMU_ERRDEF.gyro_bias_var)/reducefactor;
rng('default');
[sensor_out, sensor_err]= AddIMUErr_v002(true_inp, ErrDefs,  StateModel, STD0, sampleInterval);
samples(:, 2:7)=sensor_out';
fileID = fopen(noisyImuFile,'w');
fimu=fopen(trueImuFile,'r');
h =['% timestamp(sec), $a_{is}^s$(m/s^2), $Omega_{is}^s$(rad/sec), error in a(m/s^2), error in \omega(rad/sec)'];
fprintf(fileID, '%s\n', h);
fprintf(fileID,'%f %f %f %f %f %f %f %f %f %f %f %f %f\n',[samples(:,1)'; sensor_out; sensor_err(1:6,:)]);
fclose(fileID);
% disp(sensor_out)
% disp(sensor_err)
end