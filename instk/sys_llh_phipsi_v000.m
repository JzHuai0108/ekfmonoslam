%System model for a geodesic frame navigation based on Llh position mechanziation. 
%att must be either Cbn or qbn.
%Position errors are modelled as delta_Llh (not in metric).

%The script contains both phi and psi formulation which can be selected
%using mode parameter:
%mode=1 --> phi formulation
%mode=2 --> psi formulation

%Note that if the position errors are replaced with the metric errors (i.e.
%delta_lon= delta_Re/(R+h)/cosL, delta_lat=-delta_Rn/(R+h) ), the resulting
%error model will be equal to the model derived in Benson (1975) for both phi and psi)
%(of course, you should alse replace the delta_llh model with the
%delta_metric model defined in Benson)

%%Also, you should always keep in mind that you can use whatever position
%%error formulation as you wish (e.g. delta_metric, delta_theta or
%%delta_llh). In this toolkit I generally prefer delta_llh as it is easier
%%for new students. Benson uses delta_metric, Savage uses delta_theta. But,
%%keep in mind that they are all the same (not theoretically, but
%%practically)

%%imutype: a parameter to be directly passed to imu_err_defs_v000
function [STM Qd]=sys_llh_phipsi(Llh, Vn, att, acc, gyro, dt, mode, imutype)
% position   (1-3)
% velocity   (4-6)
% attitude   (7-9)

if (size(att,1)==1 || size(att,2)==1) %att is qbn
    Cbn=quat2dcm_v000(att);
else %att is dcm
    Cbn=att;
end

cL=cos(Llh(1));
sL=sin(Llh(1));
[Fc, wen_n, wie_n, g]=geoparam_v001(2, [cL 0 -sL]', Llh(3), Vn);

%%%%Navigation Error parameters
%system disturbance coefs
Nnav=zeros(9,6);
Nnav(7:9,4:6)=-Cbn; %attitude
Nnav(4:6,1:3)=Cbn; %velocity

%%%System Matrix
Anav=zeros(9);
%1. Common elements for both phi and psi
%1.1 Attitude Errors
Anav(7:9,7:9)=-skew(wen_n+ wie_n);

%1.2 Velocity Errors
Anav(6,3)=-2*g*Fc(1,1);  %effect of height on gravity
Anav(4:6,4:6)=-skew(wen_n+2*wie_n);
Anav(4:6,7:9)=skew(Cbn*acc);

%1.3 Position errors
Anav(1,3)=-Vn(1)*Fc(2,2)^2; 
Anav(1,4)=Fc(2,2);

Anav(2,1)=-Vn(2)*Fc(3,1)/cL;
Anav(2,3)=-Vn(2)*Fc(1,1)^2/cL;
Anav(2,5)=Fc(1,1)/cL;

Anav(3,6)=-1;

%2. PSI formulation
if (mode==2)
    %2.1 Velocity errors
    %effect of theta on gravity (we know the gravity in "n" frame. This
    %error comes from Cnc which is required to transform n-frame to the c-frame)
    Anav(4:5,1:2)=g*[-1 0;0 cL];
    
    %2.2 Position errors
    %Note that dV_n=dV_c+S(Vc)*theta
    Anav(1:3,1:2)=Anav(1:3,1:2)+skew(Vn)*[0 cL;-1 0;0 -sL];
end

%3. PHI formulation
%Here we have to add the transport rate error on attitude and velocity
%states. (In psi formaulation we do not have these components because they
%are known quantities (in computed frame) by definition)
if (mode==1)
    %3.1 Attitude errors
    %3.1.1 :effect of dw_ie_n
    Anav(7:9,1)=Anav(7:9,1)+[wie_n(3);0;-wie_n(1)];
    
    %3.1.2 :effect of dw_en_n
    Anav(7:9,3:5)=Anav(7:9,3:5)+[-Vn(2)*Fc(1,1)^2 0 Fc(1,1);Vn(1)*Fc(2,2)^2 -Fc(2,2) 0;Vn(2)*Fc(1,1)^2*sL/cL 0 Fc(3,1)];
    Anav(9,1)=Anav(9,1)-(Vn(2)*Fc(1,1)/cL/cL);
        
    %3.2 Velocity Errors
    %3.2.1 :effect of 2*dw_ie_n
    Anav(4:6,1)=Anav(4:6,1)+2*skew(Vn)*[wie_n(3);0;-wie_n(1)];
    %3.2.2 :effect of dw_en_n
    Anav(4:6,3:5)=Anav(4:6,3:5)+skew(Vn)*[-Vn(2)*Fc(1,1)^2 0 Fc(1,1);Vn(1)*Fc(2,2)^2 -Fc(2,2) 0;Vn(2)*Fc(1,1)^2*sL/cL 0 Fc(3,1)];
    Anav(4:6,1)=Anav(4:6,1)+skew(Vn)*[0;0;-(Vn(2)*Fc(1,1)/cL/cL)];
end


%%%%Imu error model parameters
[Aimu_d, Qimu_d, Cimu, Rimu]=imu_err_model(acc, gyro, dt, imutype);

%%%%Combine and discretize nav and imu models
Anav_d=eye(9)+dt*Anav;  %Use 1st order taylor series to discretize Anav
Qnav=Nnav*Rimu*Nnav';
Qnav_d=dt/2*(Anav_d*Qnav+Qnav*Anav_d');      %Use trapezoidal rule to discretize Rimu

STM=zeros(21,21);
STM(1:9,1:9)=Anav_d;
STM(1:9,10:end)=Nnav*Cimu*dt;
STM(10:end,10:end)=Aimu_d;

Qd=zeros(21);
Qd(1:9,1:9)=Qnav_d;
Qd(10:end,10:end)=Qimu_d;
Qd(1:9,10:end)=Nnav*Cimu*Qimu_d*dt/2;
Qd(10:end,1:9)=Qd(1:9,10:end)';

end


function [Ad, Qd, C, R]=imu_err_model(acc, gyro, dt, imutype)
%%imu output model
%%x'=Ax+w, y=[acc;gyro]'+Cx+v, where <w,w>=Q, <v,v>=R
%%x(k+1)=Ad.x(k)+w(k) <w(k),w(k)>=Qd

%NOTE THAT THE RETURNING SYSTEM MODEL(AD,QD) IS IN DISCRETE TIME, WHEREAS THE 
%OBSERVATION MODEL (C,R) IS IN CONTINUOUS TIME

if (imutype==1 || imutype==2 || imutype==3)
%Note that the following lines are model definition dependent. Even if you use another imu model (for instance a model
%withtout any scale factor error), do not modify these lines. Instead, just add a new "elseif" statement.
%I know this looks wierd. However, a completely modular structure turns out to be impractible.
    %Imu error states
    %1-3:Acc bias
    %4-6:gyro bias
    %7-9:Acc scale (ppt=part per thousand not million). ppm def can sometimes cause numerical problems
    %10-12:Gyro scale

    %load continuous time parameters
    IMU_ERRDEF=imu_err_defs_v000(imutype);

    %Discrete time model of imu error states
    acc_bias_a=exp(-dt/IMU_ERRDEF.acc_bias_Tc);
    acc_bias_w=IMU_ERRDEF.acc_bias_Tc*(IMU_ERRDEF.acc_bias_Q)/2*(1-acc_bias_a^2);  %(m/s^2)^2  (*s^2)
    gyro_bias_a=exp(-dt/IMU_ERRDEF.gyro_bias_Tc);
    gyro_bias_w=IMU_ERRDEF.gyro_bias_Tc*(IMU_ERRDEF.gyro_bias_Q)/2*(1-gyro_bias_a^2); %(rad/sec)^2 (*s^2)
    acc_scale_a=1;
    acc_scale_w=IMU_ERRDEF.acc_scale_Q*dt;     %ppt^2  (*s^2)
    gyro_scale_a=1;
    gyro_scale_w=IMU_ERRDEF.acc_scale_Q*dt;    %(ppt)^2  (*s^2)

    Ad=diag([acc_bias_a, acc_bias_a, acc_bias_a, gyro_bias_a, gyro_bias_a, gyro_bias_a, acc_scale_a, acc_scale_a, acc_scale_a, gyro_scale_a, gyro_scale_a, gyro_scale_a]);
    Qd=diag([acc_bias_w, acc_bias_w, acc_bias_w, gyro_bias_w, gyro_bias_w, gyro_bias_w, acc_scale_w, acc_scale_w, acc_scale_w, gyro_scale_w, gyro_scale_w, gyro_scale_w]);

    %Continuous time output model for the imu [acc;gyro]
    C=[eye(6), diag([acc(1),acc(2),acc(3),gyro(1),gyro(2),gyro(3)]/1000)];
    R=diag([IMU_ERRDEF.acc_vrw;IMU_ERRDEF.acc_vrw;IMU_ERRDEF.acc_vrw;IMU_ERRDEF.gyro_arw;IMU_ERRDEF.gyro_arw;IMU_ERRDEF.gyro_arw]);
else
    disp('This is an undefined error model.');
end
end
