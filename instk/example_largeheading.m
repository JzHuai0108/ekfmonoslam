%Note: although this example uses wander-azimuth concept, you should not consider this
%as a suitable wander-azimuth implementaion. This purpose of this example is to illustrate
%how large heading errors can be handled. There is an another example in the toolkit
%presenting a suitable wander-azimuth implementation.

clear;
fclose all;
%PATH='/home/yigiter/Desktop/Dropbox/INS/LargeHeading/mfiles/data/';
PATH='C:\Users\yigiter\Desktop\Dropbox\INS\LargeHeading\mfiles\data\';

%%Specify the imu error model and generate the path
[SenErrDef, IniErrDef, ObsErrDef, ini_pv, wander, dt]=sys_path_largeheading(PATH, 1);

%%Open the files
F_IMU=fopen([PATH 'imu.bin'],'rb');
F_TRU=fopen([PATH 'mnav.bin'],'rb');
F_IMUERR=fopen([PATH 'imuerr.bin'],'rb');
F_OBS_POS=fopen([PATH 'obs_pos.bin'],'rb');
F_RES=fopen([PATH 'res.bin'],'wb');
F_ERR=fopen([PATH 'err.bin'],'wb');

%propagation time for the discretized system
prp_per=round(0.33/dt);

%%Sensor error models
[dAimu, dBimu, Cimu, dDimu, sP0imu]=imu_modTI_v000(SenErrDef);  %errors models in discrete time
[Aimu, Qimu, Rimu]=dc2dc_v000(dAimu, dBimu*dBimu', dDimu*dDimu', dt, 1, []);  %models in continuous time
ximu=zeros(6,1); %imu error states

%%%%Initialize the INS
%PV is provided
pos_g=ini_pv(:,1);
vel_n=ini_pv(:,2); %note path generator is defined for body frame velocity
%Use accelerometer data to initialize the attitude
t=0;
in=0;
imu_avg=zeros(6,1);
imu_data=fread(F_IMU,8,'double'); %garbage
true_val = fread(F_TRU, 10, 'double');
imu_err=fread(F_IMUERR, 7, 'double');
avg_time=1;
while (in*dt)<avg_time
    imu_data=fread(F_IMU,8,'double'); 
    true_val = fread(F_TRU, 10, 'double');
    imu_err=fread(F_IMUERR, 7, 'double');
    imu_avg=imu_avg+imu_data(3:8);
    in=in+1;
end
imu_avg=imu_avg/in;
%Compute the roll/pitch
%[Cnb E]=align_wander(imu_avg(1:3), [0;0;-1]);
[Cnb E]=alingc_pl_v000(imu_avg(1:3), 0);
Cbn=Cnb';

%Use constant scale for all observations
[Rn, Re, g, sL, cL, WIE_E]=geoparam_v000(pos_g);
Sn=Rn+pos_g(3);
Se=(Re+pos_g(3))*cos(pos_g(1));
S=[Sn 0 0;0 Se 0;0 0 1];


%alignment mode (0=coarse, 1=fine)
%%"mode 2" is standard small angle (PHI-model) model. (I used it to verify the results.)
mode=0;

%Initial covariance
if (mode==0)
    nst=17;
    nst_nav=11;
    P=zeros(nst);
    P(10:11,10:11)=IniErrDef.wander_sP*IniErrDef.wander_sP';
    wander=[sin(0*pi/180);cos(0*pi/180)];
elseif (mode==1)
    nst=16;
    nst_nav=10;
    P=zeros(nst);
    P(10,10)=IniErrDef.att_sP(3,3)*IniErrDef.att_sP(3,3);
    %wander=[sin(105*pi/180);cos(105*pi/180)]; (use the wander provided by
    %the syspath)
    %P(10:11,10:11)=IniErrDef.wander_sP(1,1)+IniErrDef.wander_sP(2,2);
end
P(1:3,1:3)=S*IniErrDef.pos_sP*IniErrDef.pos_sP'*S';
P(4:6,4:6)=IniErrDef.vel_sP*IniErrDef.vel_sP';
P(7:8,7:8)=E(1:2,1:3)*Rimu(1:3,1:3)/avg_time*E(1:2,1:3)';
P(nst_nav+1:nst,nst_nav+1:nst)=sP0imu*sP0imu';
%Cross covarince between acc error and initial attitude
P(7:8,nst_nav+1:nst_nav+3)=E(1:2,1:3)*P(nst_nav+1:nst_nav+3,nst_nav+1:nst_nav+3);
P(nst_nav+1:nst_nav+3,7:8)=P(7:8,nst_nav+1:nst_nav+3)';



%accumulators and control flags
prp_coun=0;
STM_upd=eye(nst);
I=eye(nst);
imu_avg=zeros(6,1);

%%read the input data from the files
imu_data=fread(F_IMU,8,'double');
true_val = fread(F_TRU, 10, 'double');
imu_err=fread(F_IMUERR, 7, 'double');
obs_pos=fread(F_OBS_POS, 4, 'double');

%%Start iterations
while (~feof(F_IMU))
        imu=imu_data(3:8)-ximu;
        pr_coun=imu_data(1);
        [Cbn, vel_n, pos_g, wander]=strapdown_wander_dcm_v000(Cbn, vel_n, pos_g, wander, imu(1:3), imu(4:6), dt);
        imu_avg=imu_avg+imu;        
        prp_coun=prp_coun+1;
        
        if ((~isempty(obs_pos) && pr_coun==obs_pos(1))|| (mod(prp_coun,prp_per)==0)) %%(obs and/or covariance propagation)
            %%First propagate the covariance
            imu_avg=imu_avg/prp_coun;
            
            %%Form the discrete time navigation model
            if (mode==0 || mode==1)
                [STM Q]=sys_wander_largeheading_v000(pos_g, vel_n, Cbn, wander, imu_avg(1:3), prp_coun*dt, mode, Aimu, Qimu, Cimu, Rimu);
            elseif (mode==2)
                [STM Q]=sys_wander_dcm_v000(pos_g, vel_n, Cbn, wander, imu_avg(1:3), prp_coun*dt, Aimu, Qimu, Cimu, Rimu);
            end
            
            
            %Propagate
            P=STM*P*STM'+Q;
            prp_coun=0;
            imu_avg=zeros(6,1);
            
            %update IMU error states
            ximu=STM(nst_nav+1:nst,nst_nav+1:nst)*ximu;
            
            if (~isempty(obs_pos)) && (obs_pos(1)==pr_coun) %We have an observation to be processed
                
                if (mode==0 || mode==1) %switch from coarse/fine to small angle. The new wander angle will be fixed to the last wander estimate.
                    if (mode==0)
                        vda=P(10,10)+P(11,11);
                    elseif (mode==1)
                        vda=P(10,10);
                    end
                    if (vda<(0*pi/180)^2)
                        Cgn=[wander(2) wander(1) 0;-wander(1) wander(2) 0;0 0 1];
                        vel_g=Cgn'*vel_n;
                        mx_a=eye(17);
                        mx_a(10:11,:)=[];
                        mx_a(4:6,10:11)=-[vel_g(2) vel_g(1);-vel_g(1) vel_g(2);0 0];
                        mx_a(9,10:11)=[-wander(2) -wander(1)];
                        
                        if (mode==1)
                            mx_a(:,10)=wander(2)*mx_a(:,10)-wander(1)*mx_a(:,11);
                            mx_a(:,11)=[];
                        end
                        
                        P=mx_a*P*mx_a';
                        
                        mode=2;
                        nst_nav=9;
                        nst=15;
                        I=eye(nst);
                    end
                end
                
                
                H=zeros(3,nst); %position update
                H(:,1:3)=eye(3);
                R=S*ObsErrDef.sR*ObsErrDef.sR'*S';
                inno=S*(pos_g-obs_pos(2:4));  %position obs in m
                
                Rinno_inv=inv(H*P*H'+R);
                K=P*H'*Rinno_inv;         
                
                %correct current states
                dx=K*inno;
                %K(10:11,:)=0;
                %K(12:17,:)=0;
                pos_g=pos_g-S\dx(1:3);
                vel_n=vel_n-dx(4:6);
                mx_a=euler2dcm_v000(dx(7:9));
                Cbn=mx_a*Cbn;
                if (mode==0)    %coarse alignment
                    wander=wander-dx(10:11);
                    %%normalize wander
                    wander=wander/norm(wander);
                elseif (mode==1) %fine alignment
                    wander=[-wander(2) wander(1);wander(1) wander(2)]*[sin(dx(10));cos(dx(10))];
                end
                ximu=ximu+dx(nst_nav+1:nst);
                %update cov
                P=(I-K*H)*P*(I-K*H)'+K*R*K';

                %prepare for the next cycle
                obs_pos=fread(F_OBS_POS, 4, 'double');
                
                %Check to switch mode
                if (mode==0)    %switch from coarse to fine
                    vda=P(10,10)+P(11,11);
                    if (vda<(0*pi/180)^2)
                        mode=1;
                        nst_nav=10;
                        nst=16;
                        I=eye(nst);
                        
                        mx_a=eye(17,17);
                        mx_a(11,:)=[];
                        mx_a(10,10)=wander(2);
                        mx_a(10,11)=-wander(1);
                        P=mx_a*P*mx_a';
                    end
                end
                
            end
            
            %%Write the results (Cbg, Vb, pos_g)
            vr_a=diag(P).^0.5;
            
            imu_err=imu_err(2:7)-ximu;
            imu_std=vr_a(nst_nav+1:end);
            
            pos_err=S*(pos_g-true_val(2:4));
            pos_std=vr_a(1:3);
            
            vel_err=Cbn'*vel_n-true_val(5:7);
            mx_a=[Cbn' -Cbn'*skew(vel_n)];
            mx_b=P(4:9,4:9);
            mx_c=mx_a*mx_b*mx_a';
            vel_std=diag(mx_c).^0.5;
%             Cgn=[wander(2) wander(1) 0;-wander(1) wander(2) 0;0 0 1];
%             Cbg_true=euler2dcm_v000(true_val(8:10));
%             vel_err=Cgn'*vel_n-Cbg_true*true_val(5:7);
            
            Cbg_true=euler2dcm_v000(true_val(8:10));
            Cgn=[wander(2) wander(1) 0;-wander(1) wander(2) 0;0 0 1];
            Cbg=Cgn'*Cbn;
            mx_a=eye(3)-Cbg*Cbg_true';
            att_err=[-mx_a(2,3);mx_a(1,3);-mx_a(1,2)]; %-dcm2euler_v000(Cgn'*Cbn*Cbn_true');
            if (mode==0)
                mx_a=[Cgn' [0 0;0 0;-wander(2) wander(1)]];
            elseif (mode==1)
                mx_a=[Cgn' -[0;0;1]];
            elseif (mode==2)
                mx_a=Cgn';
            end
            mx_b=P(7:nst_nav,7:nst_nav);
            mx_c=mx_a*mx_b*mx_a';
            att_std=diag(mx_c).^0.5;
            
%             %%%Write the results (Cbn, Vn, pos_g) (only for mode 0 and
%             %%%mode1 - revise for mode2)
%             vr_a=diag(P).^0.5;
%             
%             imu_err=imu_err(2:7)-ximu;
%             imu_std=vr_a(nst_nav+1:end);
%             
%             pos_err=S*(pos_g-true_val(2:4));
%             pos_std=vr_a(1:3);
%             
%             Cbg_true=euler2dcm_v000(true_val(8:10));
%             Cgn=[wander(2) wander(1) 0;-wander(1) wander(2) 0;0 0 1];
%             Cbn_true=Cgn*Cbg_true;
%             vel_err=vel_n-Cbn_true*true_val(5:7);
%             if (mode==0)
%                 vr_a=Cbg_true*true_val(5:7);
%                 mx_a=[eye(3);[vr_a(2) vr_a(1);-vr_a(1) vr_a(2);0 0]];
%                 mx_b=P([4 5 6 10 11],[4 5 6 10 11]);
%             elseif (mode==1)
%                 vr_a=Cbn_true*true_val(5:7);
%                 mx_a=[eye(3);[vr_a(2);-vr_a(1);0]];
%                 mx_b=P([4 5 6 10],[4 5 6 10]);
%             end
%             mx_c=mx_a*mx_b*mx_a';
%             vel_std=diag(mx_c).^0.5;
%             
% 
%             mx_a=eye(3)-Cbn*Cbn_true';
%             att_err=[-mx_a(2,3);mx_a(1,3);-mx_a(1,2)]; %-dcm2euler_v000(Cgn'*Cbn*Cbn_true');
%             if (mode==0)
%                 mx_a=[eye(3) [0 0;0 0;-wander(2) -wander(1)]];
%             elseif (mode==1)
%                 mx_a=[eye(3) -[0;0;1]];
%             end
%             mx_b=P(7:nst_nav,7:nst_nav);
%             mx_c=mx_a*mx_b*mx_a';
%             att_std=diag(mx_c).^0.5;
                
            
            
            fwrite(F_ERR,[pr_coun;pos_err;vel_err;att_err;imu_err;pos_std;vel_std;att_std;imu_std],'double');            
            fwrite(F_RES, [pr_coun;pos_g; Cbn'*vel_n; dcm2euler_v000(Cgn'*Cbn);atan2(wander(1),wander(2));ximu],'double');
        end
        
        %Read the next values
        true_val = fread(F_TRU, 10, 'double');
        imu_err=fread(F_IMUERR, 7, 'double');
        imu_data=fread(F_IMU,8,'double');
end


fclose all;
err=readbin_v000([PATH 'err.bin'],31);
%imu_err=readbin_v000([PATH '\imuerr.bin'],7);
res=readbin_v000([PATH 'res.bin'],17);
tru=readbin_v000([PATH 'mnav.bin'],10);

% return;
% figure;
% plot(res(2,:), res(3,:))
% hold on
% plot(tru(2,:),(tru(3,:)),'r--');
% grid;

figure;
i=9;
plot(abs(err(1+i,:))*180/pi);
hold on;
plot(err(16+i,:)*180/pi,'r');
grid;

% figure;
% i=4;
% plot(abs(err(1+i,:)));
% hold on;
% plot(err(16+i,:),'r');
% grid;