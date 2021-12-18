%imu model parameter should be defined in cont. time
function [STM Qd]=mdl_ecef_dcm(Cbe, acc, Aimu, Qimu, Cimu, Rimu, dt)

%continuous model
[Anav N]=sys_ecef_dcm_v000(Cbe, acc, 0, []);
nst_imu=size(Aimu,1);
nst=9+nst_imu;
Ac=[Anav N*Cimu;zeros(nst_imu,9) Aimu];
Qc=[N*Rimu*N' zeros(9,nst_imu);zeros(nst_imu,9) Qimu];


%Discretize the model
mx_a = dt*[-Ac,Qc;zeros(nst),Ac'];
mx_b = expm(mx_a);
STM = mx_b(nst+1:2*nst,nst+1:2*nst)';
Qd = STM*mx_b(1:nst,nst+1:2*nst);