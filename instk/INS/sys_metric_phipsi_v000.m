%System model for geodetic navigation frame mechanization for both phi and psi parametrization.

%This script can be used for both quaternion and DCM attitude parametrizations.

%Position errors are represented as delta_r^n (position error with respect to
%earth defined in navigation frame) in meters rather than
%Delta-Llh (llh perturbation model see:sys_llh_phipsi) or Delta-Theta (Cen perturbation: to be added later)

%Delta_r^n, delta_llh and Delta_theta are assumed to be convertible to each
%other. That is why this system model can be used for any of these position
%mechanizations.

%There is no difference between the psi form of this script and the "sys_wander_psi". This is because
%theta_z does not appear anywhere in the error propagation model for the
%psi model. (Ccn*dv^c is assumed to be equal to dv^c in position error model
%and theta*g is independent of theta_z). Therefore, you can use
%sys_wander_psi instead of this script. (Also note that sys_wander_psi is a
%nonsingular implementation, whereas this script is singular at the pole as geoparam is called with arg 2)

%Returns F and Q of discrete time model. The IMU error model must be
%specified here. (Note that in the previous implementations I had seperated
%nav and imu models in different m-files. I no longer follow that approach)

%mode==1 --> PHI model
%mode==2 --> PSI model

%imutype is directly transferred to imu_err_defs_v000

function [STM Qd]=sys_metric_phipsi(Cen, h, Vn, att, acc, gyro, dt, mode, imutype, modelNo)

%%%Navigation error model
%1-3:Position error in meters
%4-6:vel errors
%7-9:Attitude errors in psi form

%N[da:dw]':Effect of IMU errors on nav states
if(nargin<10)
    modelNo=3; % random walk bias and scale factor
end
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
    cL=Cen(3,1);    %Note that I assume wander=0 for Cen. (This is not a non-singular implementation)
    %%position errors
    Anav(1:3,1:2)=Anav(1:3,1:2)-skew(Vn)*[0 Fc(1,1);-Fc(2,2) 0; 0 Fc(3,1)];
    
    %%Velocity errors
    %Effect of 2*d_wie_n
    Anav(4:6,1)=Anav(4:6,1)+2*skew(Vn)*[wie_n(3);0;-wie_n(1)]*Fc(2,2);
    %Effect of d_w_en_n
    Anav(4:6,3:5)=Anav(4:6,3:5)+skew(Vn)*[Vn(2)*Fc(1,1)^2 0 Fc(1,1);-Vn(1)*Fc(2,2)^2 -Fc(2,2) 0;-Vn(2)*Fc(1,1)^2*sL/cL 0 Fc(3,1)];
    Anav(4:6,1)=Anav(4:6,1)+skew(Vn)*[0;0;-(Vn(2)*Fc(1,1)/cL/cL)]*Fc(2,2);
    
    %%Attitude errors
    %effect of dw_ie_n
    Anav(7:9,1)=Anav(7:9,1)+[wie_n(3);0;-wie_n(1)]*Fc(2,2);
    
    %effect of dw_en_n
    Anav(7:9,3:5)=Anav(7:9,3:5)+[-Vn(2)*Fc(1,1)^2 0 Fc(1,1);Vn(1)*Fc(2,2)^2 -Fc(2,2) 0;Vn(2)*Fc(1,1)^2*sL/cL 0 Fc(3,1)];
    Anav(9,1)=Anav(9,1)-(Vn(2)*Fc(1,1)/cL/cL)*Fc(2,2);
end

%PSI model
if (mode==2)
    %Velocity errors
    Anav(4:5,1:2)=Anav(4:5,1:2)+g*[-Fc(2,2) 0;0 -Fc(1,1)];    %effect of horizontal position error on gravity
end

%%%%Imu error model parameters
[Aimu_d, Qimu_d, Cimu, Rimu]=imu_err_model_v001(acc, gyro, dt, imutype, modelNo);

%%Convert continuous time Anav into discrete time
%I used van loan's method here to perform conversion so that I can obtain a
%correct Qnav_d. If you use a 1st order taylor series approximation, do not
%forget to correct the cross-correlation effects caused by arw/vrw. (In,
%the final implementation, it would be better to manually compute each
%element rather than using this method)
mx_a=dt*[-Anav,Nnav*Rimu*Nnav';zeros(9),Anav'];
mx_b = expm(mx_a);
Anav_d = mx_b(10:18,10:18)';
Qnav_d = Anav_d*mx_b(1:9,10:18);

%%Combine everything
nst_imu=size(Aimu_d);
STM=zeros(9+nst_imu);
Qd=zeros(9+nst_imu);

STM(1:9,1:9)=Anav_d;
STM(1:9,10:end)=Nnav*Cimu*dt;
STM(10:end,10:end)=Aimu_d;

Qd(1:9,1:9)=Qnav_d;
Qd(10:end,10:end)=Qimu_d;
Qd(1:9,10:end)=Nnav*Cimu*Qimu_d*dt/2;   %not necessary, but let's keep this
Qd(10:end,1:9)=Qd(1:9,10:end)';
end