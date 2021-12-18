%imu model parameter should be defined in cont. time
function [STM Qd]=mdl_ned_dcm(pos_n, vel_n, Cbn, acc, Aimu, Qimu, Cimu, Rimu, dt)

%continious model
[Anav N]=sys_ned_dcm_v000(pos_n, vel_n, Cbn, acc, 0, []);
nst_imu=size(Aimu,1);
nst=9+nst_imu;
Ac=[Anav N*Cimu;zeros(nst_imu,9) Aimu];
Qc=[N*Rimu*N' zeros(9,nst_imu);zeros(nst_imu,9) Qimu];


%Discretize the model
mx_a = dt*[-Ac,Qc;zeros(nst),Ac'];
mx_b = expm(mx_a);
STM = mx_b(nst+1:2*nst,nst+1:2*nst)';
Qd = STM*mx_b(1:nst,nst+1:2*nst);