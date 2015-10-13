%In a C implementation these parameters would be in a header file. Matlab
%does not have preprocessor. Therefore, I use such a dummy function to define
%error model parameters.

%Note that the return values of this function is not standard at all.
%Whenever, you want to use a new model, just use a new imutype and define
%whatever you want.
%Don't delete the existing ones as there may be some example scripts using
%the existing values.

function IMU_ERRDEF=imu_err_defs(imutype)

switch (imutype)
    case 1  %Error model used in OPP
        %Accel
        IMU_ERRDEF.acc_bias_Tc=3600;       %sec
        IMU_ERRDEF.acc_bias_Q=5e-9;        %(m/s^2)^2
        IMU_ERRDEF.acc_bias_var=(3e-3)^2;  %(m/s^2)^2 : SS value=(Tc*Q/2)

        IMU_ERRDEF.acc_scale_Q=1.38e-6;    %ppt^2
        IMU_ERRDEF.acc_scale_var=0.1^2;    %ppt^2 (initial var of rw)

        IMU_ERRDEF.acc_vrw=0.00066^2;      %(m/s/sqrt(s))^2

        %Gyro
        IMU_ERRDEF.gyro_bias_Tc=3600;      %sec
        IMU_ERRDEF.gyro_bias_Q=1.3e-14;    %(rad/sec)^2
        IMU_ERRDEF.gyro_bias_var=(4.83e-6)^2;  %(rad/sec)^2 : SS value=(Tc*Q/2)

        IMU_ERRDEF.gyro_scale_Q=1.38e-6;    %ppt^2
        IMU_ERRDEF.gyro_scale_var=0.1^2;    %ppt^3 (initial var of rw)

        IMU_ERRDEF.gyro_arw=(3.63e-5)^2;   %(rad/sqrt(s))^2
    case 2 %Error model for SiIMU02 (see http://www.atlanticinertial.com/uploads/pdfs/4232-SiIMU02low.pdf)
        %Accel
        IMU_ERRDEF.acc_bias_Tc=1800;       %sec
        IMU_ERRDEF.acc_bias_Q=(0.5*10*1e-3)^2*(2/IMU_ERRDEF.acc_bias_Tc);     %0.5mg SS --> (m/s^2)^2 (compute this based on SS. value and made-up time constant)
        IMU_ERRDEF.acc_bias_var=(10*10*1e-3)^2;  %10mg --> (m/s^2)^2 (repeatability)

        IMU_ERRDEF.acc_scale_var=1.5^2;    %1500ppm --> ppt^2
        IMU_ERRDEF.acc_scale_Q=1e-6;    %ppt^2

        IMU_ERRDEF.acc_vrw=(0.5/60)^2;      %0.5m/s/sqrt(hr) --> (m/s/sqrt(s))^2

        %Gyro
        IMU_ERRDEF.gyro_bias_Tc=1800;      %sec
        IMU_ERRDEF.gyro_bias_Q=(4*pi/180/3600)^2*(2/IMU_ERRDEF.gyro_bias_Tc);    %1.5d/hr - 6.5d/hr SS --> (rad/sec)^2 (compute this based on SS. value and made-up time constant)
        IMU_ERRDEF.gyro_bias_var=(75*pi/180/3600)^2;  %50d/hr - 100d/hr --> (rad/sec)^2

        IMU_ERRDEF.gyro_scale_var=0.5^2;    %500ppm --> ppt^2 
        IMU_ERRDEF.gyro_scale_Q=1e-6;    %ppt^2

        IMU_ERRDEF.gyro_arw=(0.3*pi/180/60)^2; %0.1d/sqr(hr) - 0.5d/sqrt(hr) --> (rad/sqrt(s))^2 
    case 3 %Error model for LN200
        %%Accel
        IMU_ERRDEF.acc_bias_Tc=20*60;       %sec
        IMU_ERRDEF.acc_bias_Q=(0.05*10*1e-3)^2*(2/IMU_ERRDEF.acc_bias_Tc);     %50ug (SS) --> (m/s^2)^2 
        IMU_ERRDEF.acc_bias_var=(0.6*10*1e-3)^2;  %300ug - 1mg --> (m/s^2)^2

        IMU_ERRDEF.acc_scale_var=0.6^2;    %300ppm - 5000ppm --> ppt^2
        IMU_ERRDEF.acc_scale_Q=1e-6;    %ppt^2

        IMU_ERRDEF.acc_vrw=(50*10*1e-6)^2;      %50ug/sqrt(hz) --> (m/s/sqrt(s))^2

        %%Gyro
        IMU_ERRDEF.gyro_bias_Tc=20*60;      %sec
        IMU_ERRDEF.gyro_bias_Q=(0.35*pi/180/3600)^2*(2/IMU_ERRDEF.gyro_bias_Tc);    %0.35d/hr SS --> (rad/sec)^2
        IMU_ERRDEF.gyro_bias_var=(2*pi/180/3600)^2;  %1d/hr - 3d/hr --> (rad/sec)^2

        IMU_ERRDEF.gyro_scale_var=0.3^2;    %100ppm-500ppm --> ppt^2 
        IMU_ERRDEF.gyro_scale_Q=1e-6;    %ppt^2

        IMU_ERRDEF.gyro_arw=(0.1*pi/180/60)^2; %0.07d/sqr(hr) - 0.15d/sqrt(hr) --> (rad/sqrt(s))^2 
     case 4 %Error model for H764G
        %%Accel  
        IMU_ERRDEF.acc_bias_Tc=20*60;       %sec
        IMU_ERRDEF.acc_bias_Q=(25*10*1e-6)^2*(2/IMU_ERRDEF.acc_bias_Tc);     %50ug (SS) --> (m/s^2)^2 
        IMU_ERRDEF.acc_bias_var=(25*10*1e-6)^2;  %25ug - 1mg --> (m/s^2)^2
        IMU_ERRDEF.initacc_bias_err=25*10*1e-6;
        
        IMU_ERRDEF.acc_scale_var=(100e-3)^2;    %100ppm --> ppt^2, note another 1e-3 was accounted for in IMU error compensation   
        IMU_ERRDEF.acc_scale_Q=1e-6;    %ppt^2
        IMU_ERRDEF.acc_vrw=(8.3e-6)^2;      %8.3ug/sqrt(100hz) --> (m/s/sqrt(s))^2

        %%Gyro   
        IMU_ERRDEF.gyro_bias_Tc=20*60;      %sec
        IMU_ERRDEF.gyro_bias_Q=(0.0035*pi/180/3600)^2*(2/IMU_ERRDEF.gyro_bias_Tc);    %0.35d/hr SS --> (rad/sec)^2
        IMU_ERRDEF.gyro_bias_var=(.0035*pi/180/3600)^2;  %.0035d/hr--> (rad/sec)^2
        IMU_ERRDEF.initgyro_bias_err=.0035*pi/180/3600;% for high end, no turn on bias
        
        IMU_ERRDEF.gyro_scale_var=(5e-3)^2;    %5ppm --> ppt^2, note another 1e-3 was accounted for in IMU error compensation    
        IMU_ERRDEF.gyro_scale_Q=1e-6;    %ppt^2
        IMU_ERRDEF.gyro_arw=(0.0035*pi/180/60)^2; %0.0035d/sqr(hr) --> (rad/sqrt(s))^2 
    case 5% mems microstrain 3dm gx 3 35
        IMU_ERRDEF.acc_bias_Tc=1800;       %sec
        IMU_ERRDEF.acc_bias_Q=(.04*1.0e-2)^2*(2/IMU_ERRDEF.acc_bias_Tc); 
        IMU_ERRDEF.acc_bias_var=(.04*1.0e-2)^2; % this parameter's meaning is different from its original meaning, turn on bias variance,
        % now it means bias drift variability
        IMU_ERRDEF.initacc_bias_err=0.002*9.8;% 2mg  turn on bias
        
        IMU_ERRDEF.gyro_bias_Tc=1800;      %sec
        IMU_ERRDEF.gyro_bias_Q=(18*(pi/180.0)/3600)^2*(2/IMU_ERRDEF.gyro_bias_Tc); % 18deg/hr  
        IMU_ERRDEF.initgyro_bias_err=0.25*pi/180; %.25 deg turn on bias
        IMU_ERRDEF.gyro_bias_var=(18*(pi/180.0)/3600)^2;
   
        IMU_ERRDEF.acc_scale_var=(500.0*1.0e-3)^2; % 500 ppm*1000
        IMU_ERRDEF.acc_scale_Q=100*1e-6;    %ppt^2, Tc~3600sec
        IMU_ERRDEF.gyro_scale_var = (500*1.0e-3)^2; % 500 ppm*1000
        IMU_ERRDEF.gyro_scale_Q=100*1e-6;    %ppt^2, Tc~3600sec

        IMU_ERRDEF.acc_vrw= (80*1.0e-5)^2; % 80 ug with 1 Hz    
        IMU_ERRDEF.gyro_arw = (0.03*pi/180)^2; % 0.03 deg/sec/sqrt(Hz) 
    case 6% mems steval mki062v2, assumed to be 5 times worse than 3dm gx3-35
        IMU_ERRDEF.acc_bias_Tc=1800;       %sec
        IMU_ERRDEF.acc_bias_Q=(5*0.04*1.0e-2)^2*(2/IMU_ERRDEF.acc_bias_Tc); 
        IMU_ERRDEF.acc_bias_var=(5*.04*1.0e-2)^2; % this parameter's meaning is different from its original meaning, turn on bias variance,
        % now it means bias drift variability
        IMU_ERRDEF.initacc_bias_err=5*0.002*9.8;% 5x2mg  turn on bias
        
        IMU_ERRDEF.gyro_bias_Tc=1800;      %sec
        IMU_ERRDEF.gyro_bias_Q=(5*18*(pi/180.0)/3600)^2*(2/IMU_ERRDEF.gyro_bias_Tc); % 5x18deg/hr  
        IMU_ERRDEF.initgyro_bias_err=5*0.25*pi/180; %5x.25 deg turn on bias
        IMU_ERRDEF.gyro_bias_var=(5*18*(pi/180.0)/3600)^2;
   
        IMU_ERRDEF.acc_scale_var=(500.0*1.0e-3)^2; % 500 ppm*1000
        IMU_ERRDEF.acc_scale_Q=25*100*1e-6;    %ppt^2, Tc~3600sec
        IMU_ERRDEF.gyro_scale_var = (500*1.0e-3)^2; % 500 ppm*1000
        IMU_ERRDEF.gyro_scale_Q=25*100*1e-6;    %ppt^2, Tc~3600sec

        IMU_ERRDEF.acc_vrw= (5*80*1.0e-5)^2; % 5*80 ug with 1 Hz    
        IMU_ERRDEF.gyro_arw = (5*0.03*pi/180)^2; % 5*0.03 deg/sec/sqrt(Hz) 
    case 7 % EPSON m-g362pdc1 mems imu
        IMU_ERRDEF.acc_bias_Tc=1800;       %sec
        IMU_ERRDEF.acc_bias_Q=(1.0e-3)^2*(2/IMU_ERRDEF.acc_bias_Tc);         
        IMU_ERRDEF.acc_bias_var=(1e-3)^2; % this parameter's means bias drift variability
        IMU_ERRDEF.initacc_bias_err=8e-2;% 8mg  turn on bias
        
        IMU_ERRDEF.gyro_bias_Tc=1800;      %sec
        IMU_ERRDEF.gyro_bias_Q=(3*(pi/180.0)/3600)^2*(2/IMU_ERRDEF.gyro_bias_Tc); % 3 deg/hr 
        IMU_ERRDEF.initgyro_bias_err=0.5*pi/180; %.5 deg turn on bias
        IMU_ERRDEF.gyro_bias_var=(3*(pi/180.0)/3600)^2; %3 deg/hr
   
        IMU_ERRDEF.acc_vrw= (0.04/60)^2; % 0.04 m/s/sqrt(hr)    
        IMU_ERRDEF.gyro_arw = (10*0.1*pi/180/60)^2; % 0.1 deg/sqrt(hr), but in tests I found it is much larger, so multiply it by 10
end

