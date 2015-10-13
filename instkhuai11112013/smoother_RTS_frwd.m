%Bf formula is more flexible than the original rts formula. (Essentially,
%most people implements bf while they try to implement rts without realizing
%this fact!)

clear;
fclose all;
randn('state',[112233;445566]); %make non-random (required for debugging)
PATH='C:\Work\Benisil';

%%Specify the imu error model and generate the path
[SenErrDef, IniErrDef, ObsErrDef, ini_pva, dt]=sys_path(PATH);
nst=15;

%%Open the files
%Filter
F_IMU1=fopen([PATH '\imu.bin'],'rb');
F_TRU1=fopen([PATH '\mnav.bin'],'rb');
F_IMUERR1=fopen([PATH '\imuerr.bin'],'rb');
F_OBS_POS1=fopen([PATH '\obs_pos.bin'],'rb');
F_RES1=fopen([PATH '\res1.bin'],'wb');

%smoother
F_PROP=fopen([PATH '\smo_prop.bin'],'wb');
F_UPDT=fopen([PATH '\smo_updt.bin'],'wb');

%%Program parameters
prp_per=round(0.33/dt); %propagation time for the discretized system

%%Sensor error models
[dAimu, dBimu, Cimu, dDimu, sP0imu]=imu_modTI_v000(SenErrDef);  %errors models in discrete time
[Aimu, Qimu, Rimu]=dc2dc_v000(dAimu, dBimu*dBimu', dDimu*dDimu', dt, 1, []);  %models in continuous time
%%%%%%%%Rimu=Rimu/dt; %%correction for the notation difference. (we need a standard notation for mems errors!!!!!)

%Initial covariance
P1=zeros(nst);
P1(1:3,1:3)=IniErrDef.pos_sP*IniErrDef.pos_sP';
P1(4:6,4:6)=IniErrDef.vel_sP*IniErrDef.vel_sP';
P1(7:9,7:9)=IniErrDef.att_sP*IniErrDef.att_sP';
P1(10:nst,10:nst)=sP0imu*sP0imu';

%imu error states
ximu1=zeros(6,1);

%INS registers
pos_n1=ini_pva(:,1);
att_n1=ini_pva(:,3);
Cbn1=euler2dcm_v000(att_n1);
vel_n1=Cbn1*ini_pva(:,2); %note path generator is defined for body frame velocity

%accumulators and control flags
avg_imu1=zeros(6,1);

pr_coun1=0;prp_coun1=0;
STM_upd=eye(nst);

write_flag1=0;
I=eye(nst);

%%read the input data from the files
imu_data1=fread(F_IMU1,8,'double'); %garbage
imu_data1=fread(F_IMU1,8,'double');
true_val1 = fread(F_TRU1, 10, 'double');
imu_err1=fread(F_IMUERR1, 7, 'double');

obs_pos1=fread(F_OBS_POS1, 4, 'double');

%%Start iterations
while (~feof(F_IMU1))
        imu1=imu_data1(3:8)-ximu1;
        [Cbn1, vel_n1, pos_n1]=strapdown_ned_dcm_v000(Cbn1, vel_n1, pos_n1, imu1(1:3), imu1(4:6), dt);
        avg_imu1=avg_imu1+imu1;
        
        pr_coun1=pr_coun1+1;
        prp_coun1=prp_coun1+1;
        
        if (mod(prp_coun1,prp_per)==0 && prp_coun1~=0) %propagate the covariance and estimates
            avg_imu1=avg_imu1/prp_coun1;
            
            %%Form the discrete time navigation model
            [STM1 Q1]=mdl_ned_dcm_v000(pos_n1, vel_n1, Cbn1, avg_imu1(1:3), Aimu, Qimu, Cimu, Rimu, prp_coun1*dt);
            
            %Propagate
            P1=STM1*P1*STM1'+Q1;
            STM_upd=STM1*STM_upd;
            prp_coun1=0;
            avg_imu1=zeros(6,1);
            
            %update IMU error states
            ximu1=STM1(10:nst,10:nst)*ximu1;
            
            write_flag1=1;
            
            %smoother stuff;
            smo_STM=STM_upd(:);
            smo_p_pre=P1(:);
            STM_upd=eye(nst);
        end
        
        if (~isempty(obs_pos1) && pr_coun1==obs_pos1(1))
            %first check whether the covariance is uptodate
            if (prp_coun1~=0)
                avg_imu1=avg_imu1/prp_coun1;
                [STM1 Q1]=mdl_ned_dcm_v000(pos_n1, vel_n1, Cbn1, avg_imu1(1:3), Aimu, Qimu, Cimu, Rimu, prp_coun1*dt);
                P1=STM1*P1*STM1'+Q1;
                STM_upd=STM1*STM_upd;
                prp_coun1=0;
                avg_imu1=zeros(6,1);
                ximu1=STM1(10:nst,10:nst)*ximu1;

                %smoother stuff;
                smo_STM=STM_upd(:);
                smo_p_pre=P1(:);
                STM_upd=eye(15);
            end
            
            smo_dx_inp=0;
            smo_dP_inp=0;   %I do not need this.
            while (~isempty(obs_pos1) && (pr_coun1==obs_pos1(1)))   %there may be more than one observation at any instant
                H1=zeros(3,nst); %position update
                H1(:,1:3)=eye(3);
                R1=ObsErrDef(1).sR*ObsErrDef(1).sR';
                inno1=pos_n1-obs_pos1(2:4);

                Rinno_inv1=inv(H1*P1*H1'+R1);
                K1=P1*H1'*Rinno_inv1;         

                %smoother stuff
                smo_dx_inp=smo_dx_inp+STM_upd'*H1'*Rinno_inv1*inno1(:);
                smo_dP_inp=smo_dP_inp+STM_upd'*H1'*Rinno_inv1*H1*STM_upd; %I do not need this.(RTS uses P_filt-P_pre as the input)
                
                %correct current states
                dx1=K1*inno1;
                pos_n1=pos_n1-dx1(1:3);
                vel_n1=vel_n1-dx1(4:6);
                mx_a=euler2dcm_v000(dx1(7:9));
                Cbn1=mx_a*Cbn1;
                ximu1=ximu1+dx1(10:nst);
                %update cov
                P1=(I-K1*H1)*P1*(I-K1*H1)'+K1*R1*K1';

                %prepare for the next cycle
                STM_upd=(I-K1*H1)*STM_upd;
                obs_pos1=fread(F_OBS_POS1, 4, 'double');

                write_flag1=1;
            end
            %Smoother stuff
            P1_pre=reshape(smo_p_pre,15,15);
            smo_dx=P1_pre*smo_dx_inp;
            %smo_dx=smo_dx_inp;
            fwrite(F_UPDT,[pr_coun1;smo_dx],'double');
            STM_upd=eye(15);
        end
        
        %True values
        true_val1 = fread(F_TRU1, 10, 'double');
        imu_err1=fread(F_IMUERR1, 7, 'double');
        
        %%%write the results
        %read the true values
        if (write_flag1)
            err1=zeros(nst,1);
            err1(1:3)=pos_n1-true_val1(2:4);
            Cbn_true=euler2dcm_v000(true_val1(8:10));
            err1(4:6)=vel_n1-Cbn_true*true_val1(5:7);
            err1(7:9)=-dcm2euler_v000(Cbn1*Cbn_true');
            err1(10:nst)=imu_err1(2:7)-ximu1;
            
            fwrite(F_RES1,[pr_coun1;err1;diag(P1).^0.5],'double');
            write_flag1=0;
            
            %smoother stuff
            smo_p_upd=P1(:);
            fwrite(F_PROP,[pr_coun1;smo_p_pre;smo_p_upd;smo_STM],'double');
        end
        
        %read the next imu output
        imu_data1=fread(F_IMU1,8,'double');
end


fclose all;
res1=readbin_v000([PATH '\res1.bin'],31);

figure;
plot(res1(1,:)/100,abs(res1(2,:))*6e6)
hold on
plot(res1(1,:)/100,(res1(17,:))*6e6,'--');
title('lat. error');
ylabel('m');
xlabel('sec');
legend('filter error', 'filter std');
grid;
return;


figure;
yy=6;
plot(res1(1,:)/100,abs(res1(1+yy,:)))
hold on
plot(res1(1,:)/100,(res1(16+yy,:)),'--');
xlabel('sec');
legend('filter', 'fls');
grid;

figure;
yy=2;
plot(res1(1,:)/100,abs(res1(1+yy,:))*6e6)
hold on
plot(res1(1,:)/100,(res1(16+yy,:))*6e6,'--');
xlabel('sec');
legend('filter', 'fls');
grid;

