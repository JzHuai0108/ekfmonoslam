%%Examples showing the usage of coarse align methods

clear;
%Parameter Values
eul=[-25;46;143]*pi/180;
Cbn=euler2dcm_v000(eul);
Cpn=euler2dcm_v000([0;0;eul(3)]);
Cbp=euler2dcm_v000([eul(1);eul(2);0]);
Llh=[35*pi/180;60*pi/180;900];

[Rn, Re, g, sL, cL, WIE_E]=geoparam_v000(Llh);
acc=Cbn'*[0;0;-g];
gyro=Cbn'*[WIE_E*cL;0;WIE_E*-sL];

acc_err=[0.2;-0.5;-0.01];
gyro_err=[-1;2;-3]*1e-5;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%Estimation Results

%% Cpb_est (heading unknown)
[Cpb_est Ep]=alingc_pl_v000(acc+acc_err, 0);
SP=Cbp*Cpb_est-eye(3);	%estimation error
[Ep*acc_err [SP(3,2);-SP(3,1);SP(2,1)]]		%comparison of true error with linear error model

%% Cnp_est (using the correct Cbp)
[Cnp_est Ehg Ehp]=alingc_hd_v000(gyro+gyro_err, Llh(1), Cbp);
SH=Cpn*Cnp_est-eye(3);
[Ehg*gyro_err SH(2,1)]

%% Cnb_est1 (with true heading)
[Cnb_est1 Ep1]=alingc_pl_v000(acc+acc_err, eul(3));
SN1=Cbn*Cnb_est1-eye(3);
[Ep1*acc_err [SN1(3,2);-SN1(3,1);SN1(2,1)]]

%% Cnb_est (heading is estimated from gyro)
[Cnb_est E]=alingc_2s_v000(acc+acc_err, gyro+gyro_err, Llh);
SN=Cbn*Cnb_est-eye(3);
[E*[acc_err;gyro_err] [SN(3,2);-SN(3,1);SN(2,1)]]

%% Cnb_est2 (Std gyro compass method)
Cnb_est2=alingc_std_v001(acc,gyro,Llh, 1)

%% Cnb_est3 (Using ref vector)
[Cnb_est3 Ea]=alingc_std_v001(acc+acc_err,gyro+gyro_err,[], 2);
SA=Cbn*Cnb_est3-eye(3);
[Ea*[acc_err;gyro_err] [SA(3,2);-SA(3,1);SA(2,1)]]
