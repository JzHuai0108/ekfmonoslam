%2 Step coarse align method.
%if Cbn_est=(I-S(e))Cbn, then e=E[acc_err;gyro_err]
function [Cnb E]=alingc_2st(acc, gyro, Llh)

%Compute the Cpb using accelerometers
[Cpb Ep]=alingc_pl_v000(acc,0);

%Compute the heading transformation from gyro
[Cnp Ehg Ehp]=alingc_hd_v000(gyro, Llh(1), Cpb');

%Final Transformation
Cnb=Cpb*Cnp;

%Error matrix
E=[Cnp'*Ep zeros(3)]+[zeros(2,6);Ehp*Ep Ehg];
