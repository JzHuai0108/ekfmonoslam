% This function assumes bias and scale factor errors of gyros and
% accelerometers are random constant, compare with sys_metric_phipsi_v000
% where they are modeled as random walk processes. 
% Both models can be used in the %forwardfilterxxxx.m.
function [STM Qd]=sys_metric_phipsi(Cen, h, Vn, att, acc, gyro, dt, mode, imutype)

%%%Navigation error model
%1-3:Position error in meters
%4-6:vel errors
%7-9:Attitude errors in psi form

%N[da:dw]':Effect of IMU errors on nav states

Anav=zeros(9); %pos, vel, attitude
Nnav=zeros(9,6);

[Fc, wen_n, wie_n, g]=geoparam_v001(2, Cen(:,3), h, Vn);

if (size(att,1)==1 || size(att,2)==1) %att is qbn
    Cbn=quat2dcm_v000(att);
else %att is dcm
    Cbn=att;
end

%%%Navigation System model
%%Common terms for both phi and psi
%position errors
Anav(1:3,1:3)=-skew(wen_n);
Anav(1:3,4:6)=eye(3);

%Velocity Errors
Anav(6,3)=2*g*Fc(1,1);  %effect of height on gravity (note that, drz is the negative of dh. Therefore, this is positive)
Anav(4:6,4:6)=-skew(wen_n+2*wie_n);
Anav(4:6,7:9)=skew(Cbn*acc);

Nnav(4:6,1:3)=Cbn; %%effect of accelerometer errors

%Attitude errors
Anav(7:9,7:9)=-skew(wen_n+ wie_n);
Nnav(7:9,4:6)=-Cbn; %Effect of gyroscope errors

%PHI model
if (mode==1)
    sL=-Cen(3,3);
    cL=Cen(3,1);    %Note that I assume wander=0 for Cen. (This is not a singular implementation)
    %%position errors
    Anav(1:3,1:2)=Anav(1:3,1:2)-skew(Vn)*[0 Fc(1,1);-Fc(2,2) 0; 0 Fc(3,1)];
    
    %%Velocity errors
    %Effect of 2*d_wie_n
    Anav(4:6,1)=Anav(4:6,1)+2*skew(Vn)*[wie_n(3);0;-wie_n(1)]/Fc(2,2);
    %Effect of d_w_en_n
    Anav(4:6,3:5)=Anav(4:6,3:5)+skew(Vn)*[Vn(2)*Fc(1,1)^2 0 Fc(1,1);-Vn(1)*Fc(2,2)^2 -Fc(2,2) 0;-Vn(2)*Fc(1,1)^2*sL/cL 0 Fc(3,1)];
    Anav(4:6,1)=Anav(4:6,1)+skew(Vn)*[0;0;-(Vn(2)*Fc(1,1)/cL/cL)]/Fc(2,2);
    
    %%Attitude errors
    %effect of dw_ie_n
    Anav(7:9,1)=Anav(7:9,1)+[wie_n(3);0;-wie_n(1)]/Fc(2,2);
    
    %effect of dw_en_n
    Anav(7:9,3:5)=Anav(7:9,3:5)+[-Vn(2)*Fc(1,1)^2 0 Fc(1,1);Vn(1)*Fc(2,2)^2 -Fc(2,2) 0;Vn(2)*Fc(1,1)^2*sL/cL 0 Fc(3,1)];
    Anav(9,1)=Anav(9,1)-(Vn(2)*Fc(1,1)/cL/cL)/Fc(2,2);
end

%PSI model
if (mode==2)
    %Velocity errors
    Anav(4:5,1:2)=Anav(4:5,1:2)+g*[-Fc(2,2) 0;0 -Fc(1,1)];    %effect of horizontal position error on gravity
end

%%Combine everything
nst_imu=12;
STM=zeros(9+nst_imu);
Qd=zeros(9+nst_imu);
Gk=zeros(9+nst_imu,6);
STM(1:9,1:9)=Anav;
STM(4:6,10:18)=[Cbn, zeros(3), Cbn*diag(acc)/1000];
STM(7:9,13:21)=[-Cbn, zeros(3),-Cbn*diag(gyro)/1000];
STM=eye(21)+STM*dt; % FIRST ORDER TAYLOR APPROXIMATION

Gk( 4: 6, 1: 3)	= Cbn;		% noise driving accelerometer 
Gk( 7: 9, 4: 6)	=-Cbn;		% noise driving gyro  
IMU_ERRDEF=imu_err_defs_v000(imutype);
R=diag([IMU_ERRDEF.acc_vrw;IMU_ERRDEF.acc_vrw;IMU_ERRDEF.acc_vrw;IMU_ERRDEF.gyro_arw;IMU_ERRDEF.gyro_arw;IMU_ERRDEF.gyro_arw]);
Qd=Gk*R*dt*Gk';
end
