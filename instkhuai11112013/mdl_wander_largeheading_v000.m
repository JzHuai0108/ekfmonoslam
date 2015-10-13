%imu model parameter should be defined in cont. time
function [STM Qd]=mdl_wander_largeheading(pos_n, vel_n, Cbn, wander, acc, Aimu, Qimu, Cimu, Rimu, dt)

%continuous model
[Anav N]=sys_wander_largeheading_v000(pos_n, vel_n, Cbn, wander, acc, 0, []);
nst_imu=size(Aimu,1);
nst=11+nst_imu;
Ac=[Anav N*Cimu;zeros(nst_imu,11) Aimu];
Qc=[N*Rimu*N' zeros(11,nst_imu);zeros(nst_imu,11) Qimu];


%Discretize the model
mx_a = dt*[-Ac,Qc;zeros(nst),Ac'];
mx_b = expm(mx_a);
STM = mx_b(nst+1:2*nst,nst+1:2*nst)';
Qd = STM*mx_b(1:nst,nst+1:2*nst);
