fclose all;
clear;

%%%% Calibration
%g0=9.80665; %standard garvity value
g=9.80801; %real gravity in the lab


%Sync and interpolate the data
% DIRNAME='C:\Benisil\PDANav\CalibTest1\';
% INPFILES=[[DIRNAME '07-16-00-12-acclog.bin ']; [DIRNAME '07-16-00-12-gyrolog.bin']];
% interp_imu_v001(INPFILES, [DIRNAME 'imu.bin'], [5 2 3 4 5;5 2 3 4 5], 2);

imudata=readbin_v000([DIRNAME 'imu.bin'],7);
bounind=[35 95;105 160;180 250;265 315;355 440;455 510;520 585;600 690;700 790;800 890;900 960;
    970 1070;1120 1175;1190 1270;1285 1345;1355 1430;1445 1500;1515 1570;1650 1720;1740 1790;
    1870 1945;1960 2050;2065 2115;2170 2260;2280 2335;2365 2425;2435 2510;2520 2575;2585 2645;
    2655 2720;2725 2785];


val_set=[0; 1; 150*pi/180; 300*pi/180; 450*pi/180; 500*pi/180];
mx_a=[0 1 0;1 0 0; 0 0 -1];
mx_b=[0 -1 0;-1 0 0; 0 0 1];
C=diagmat_v000(mx_a,mx_b,1);
[dout din]=cp_sensorio_v000(imudata(2:end,:), bounind, val_set, C);
din(1:3,:)=din(1:3,:)*g;
%dout(1:3,:)=dout(1:3,:)*g;

in_st=find(sum(din(4:6,:),1)==0);
in_rt=find(sum(din(4:6,:),1)~=0);

%Acc calib matrix and bias
[Aacc bacc]=cp_affine_v000(dout(1:3,in_st),din(1:3,in_st));
dacc_calib=Aacc\(dout(1:3,:)-bacc*ones(1,size(din,2)));

%Gyro bias and g-dependence
[Bgyro bgyro1]=cp_affine_v000(dout(4:6,in_st),din(1:3,in_st));

%Gyro Calib matrix and bias
gdat=dout(4:6,in_rt)-(Bgyro*dacc_calib(1:3,in_rt)+bgyro1*ones(1,length(in_rt)));
[Agyro bgyro2]=cp_affine_v000(gdat,din(4:6,in_rt));

% %Gyro Calib matrix g-dependence and bias
% [mx_a vr_a]=cp_affine_v000(dout(4:6,:),din);
dgyro_calib=Agyro\(dout(4:6,:)-(Bgyro*dacc_calib(1:3,:))-(bgyro1+bgyro2)*ones(1,size(dout,2)));


return;
figure;
plot(dgyro_calib(1,:),'+');
hold on
plot(din(4,:),'r+')
plot(dout(4,:),'k+')
grid;

figure;
plot(dacc_calib(2,in_st),'+');
hold on
plot(din(2,in_st),'r+')
plot(dout(2,in_st)*g,'k+')
grid;