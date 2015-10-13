clear;
fclose all;
randn('state',666); %make non-random (required for debugging)
PATH='C:\Work\Benisil';

%%Specify the imu error model and generate the path
[SenErrDef, IniErrDef, ObsErrDef, ini_pva, dt]=sys_path(PATH);
%copy the files for the second(lag state) and third(buffer) ins implementations
copyfile([PATH '\imu.bin'], [PATH '\imu1.bin']);
copyfile([PATH '\obs_pos.bin'], [PATH '\obs_pos1.bin']);
copyfile([PATH '\mnav.bin'], [PATH '\mnav1.bin']);
copyfile([PATH '\imuerr.bin'], [PATH '\imuerr1.bin']);

%%Open the files
F_IMU1=fopen([PATH '\imu.bin'],'rb');
F_IMU2=fopen([PATH '\imu1.bin'],'rb');
F_TRU1=fopen([PATH '\mnav.bin'],'rb');
F_TRU2=fopen([PATH '\mnav1.bin'],'rb');
F_IMUERR1=fopen([PATH '\imuerr.bin'],'rb');
F_IMUERR2=fopen([PATH '\imuerr1.bin'],'rb');
F_OBS_POS1=fopen([PATH '\obs_pos.bin'],'rb');
F_OBS_POS2=fopen([PATH '\obs_pos1.bin'],'rb');
F_RES1=fopen([PATH '\res1.bin'],'wb');
F_RES2=fopen([PATH '\res2.bin'],'wb');
F_RES3=fopen([PATH '\res3.bin'],'wb');

%%Program parameters
delay=11/dt; %delay for the fixed lag 
prp_per=0.5/dt; %propagation time for the discretized system
dx_upd_per=20/dt; %accumulator correction period 

%%Sensor error models
[dAimu, dBimu, Cimu, dDimu, sP0imu]=imu_modTI_v000(SenErrDef);  %errors models in discrete time
[Aimu, Qimu, Rimu]=dc2dc_v000(dAimu, dBimu*dBimu', dDimu*dDimu', dt, 1, []);  %models in continious time
%%%%%%%%Rimu=Rimu/dt; %%correction for the notation difference. (we need a standard notation for mems errors!!!!!)

%Initial covariance
P1=zeros(15);
P1(1:3,1:3)=IniErrDef.pos_sP*IniErrDef.pos_sP';
P1(4:6,4:6)=IniErrDef.vel_sP*IniErrDef.vel_sP';
P1(7:9,7:9)=IniErrDef.att_sP*IniErrDef.att_sP';
P1(10:15,10:15)=sP0imu*sP0imu';
P2=P1;

%imu error states
ximu1=zeros(6,1);
ximu2=zeros(6,1); %buffer ins

%INS registers
pos_n1=ini_pva(:,1);
att_n1=ini_pva(:,3);
Cbn1=euler2dcm_v000(att_n1);
vel_n1=Cbn1*ini_pva(:,2); %note path generator is defined for body frame velocity

pos_n2=pos_n1; %buffer ins
vel_n2=vel_n1;
Cbn2=Cbn1; 

%accumulators and control flags
avg_imu1=zeros(6,1);
avg_imu2=zeros(6,1);

pr_coun1=0;prp_coun1=0;
pr_coun2=0;prp_coun2=0;
dx_coun=0;
dx_flag=1;
dx_last_up=-1;

dx12_accum=zeros(15,1);
dx11_accum=zeros(15,1);
dP12_accum=zeros(15);
dP11_accum=zeros(15);
STM12=eye(15);
STM11=eye(15);
STM_upd=eye(15);

write_flag1=0;
write_flag2=0;

I=eye(15);

%%read the input data from the files
imu_data1=fread(F_IMU1,8,'double'); %garbage
imu_data1=fread(F_IMU1,8,'double');
imu_data2=fread(F_IMU2,8,'double'); %garbage
imu_data2=fread(F_IMU2,8,'double');
true_val1 = fread(F_TRU1, 10, 'double');
imu_err1=fread(F_IMUERR1, 7, 'double');
true_val2 = fread(F_TRU2, 10, 'double');
imu_err2=fread(F_IMUERR2, 7, 'double');

obs_pos1=fread(F_OBS_POS1, 4, 'double');
obs_pos2=fread(F_OBS_POS2, 4, 'double');

%%Start iterations
while (~feof(F_IMU2))
    if (~feof(F_IMU1)) %ins1 computations (innovation INS)
        imu1=imu_data1(3:8)-ximu1;
        [Cbn1, vel_n1, pos_n1]=strapdown_ned_dcm_v000(Cbn1, vel_n1, pos_n1, imu1(1:3), imu1(4:6), dt);
        avg_imu1=avg_imu1+imu1;
        
        pr_coun1=pr_coun1+1;
        prp_coun1=prp_coun1+1;
        dx_coun=dx_coun+1;
        
        if (mod(prp_coun1,prp_per)==0 && prp_coun1~=0) %propagate the covariance and estimates
            avg_imu1=avg_imu1/prp_coun1;
            
            %%Form the discrete time navigation model
            [STM1 Q1]=mdl_ned_dcm_v000(pos_n1, vel_n1, Cbn1, avg_imu1(1:3), Aimu, Qimu, Cimu, Rimu, prp_coun1*dt);
            
            %Propagate
            P1=STM1*P1*STM1'+Q1;
            STM12=STM1*STM12;
            STM11=STM1*STM11;
            STM_upd=STM1*STM_upd;
            prp_coun1=0;
            avg_imu1=zeros(6,1);
            
            %update IMU error states
            ximu1=STM1(10:15,10:15)*ximu1;
            
            write_flag1=1;
        end
        
        
        while (~isempty(obs_pos1) && (pr_coun1==obs_pos1(1)))
            %first check wheter the covariance is uptodate
            if (prp_coun1~=0)
                avg_imu1=avg_imu1/prp_coun1;
                [STM1 Q1]=mdl_ned_dcm_v000(pos_n1, vel_n1, Cbn1, avg_imu1(1:3), Aimu, Qimu, Cimu, Rimu, prp_coun1*dt);
                P1=STM1*P1*STM1'+Q1;
                STM12=STM1*STM12;
                STM11=STM1*STM11;
                STM_upd=STM1*STM_upd;
                prp_coun1=0;
                avg_imu1=zeros(6,1);
                ximu1=STM1(10:15,10:15)*ximu1;
                
            end
            
            H1=zeros(3,15); %position update
            H1(:,1:3)=eye(3);
            R1=ObsErrDef(1).sR*ObsErrDef(1).sR';
            inno1=pos_n1-obs_pos1(2:4);
            
            Rinno_inv1=inv(H1*P1*H1'+R1);
            K1=P1*H1'*Rinno_inv1;         
            
            %correct current states
            dx1=K1*inno1;
            pos_n1=pos_n1-dx1(1:3);
            vel_n1=vel_n1-dx1(4:6);
            mx_a=euler2dcm_v000(dx1(7:9));
            Cbn1=mx_a*Cbn1;
            ximu1=ximu1+dx1(10:15);
            %update cov
            P1=(I-K1*H1)*P1*(I-K1*H1)'+K1*R1*K1';

            %correct lagged states
            vr_a=STM12'*H1'*Rinno_inv1*inno1;
            dx12_accum=dx12_accum+vr_a;
            vr_a=STM11'*H1'*Rinno_inv1*inno1;
            dx11_accum=dx11_accum+vr_a;
            dP12_accum=dP12_accum+STM12'*H1'*Rinno_inv1*H1*STM12;
            dP11_accum=dP11_accum+STM11'*H1'*Rinno_inv1*H1*STM11;

            %prepare for the next cycle
            STM12=(I-K1*H1)*STM12;
            STM11=(I-K1*H1)*STM11;
            STM_upd=(I-K1*H1)*STM_upd;
            obs_pos1=fread(F_OBS_POS1, 4, 'double');
            
            write_flag1=1;
            
            %%dx_accumulator reset
            if ((dx_coun>dx_upd_per) && (dx_flag))
                STM11=eye(15);
                dx11_accum=zeros(15,1);
                dP11_accum=zeros(15);
                dx_flag=0;
                dx_coun=0;
                dx_last_up=pr_coun1;
            end
        end
        
        %True values
        true_val1 = fread(F_TRU1, 10, 'double');
        imu_err1=fread(F_IMUERR1, 7, 'double');
        
        %%%write the results
        %read the true values
        if (write_flag1)
            err1=zeros(15,1);
            err1(1:3)=pos_n1-true_val1(2:4);
            err1(4:6)=vel_n1-true_val1(5:7);
            err1(7:9)=-dcm2euler_v000(Cbn1*euler2dcm_v000(true_val1(8:10))');
            err1(10:15)=imu_err1(2:7)-ximu1;
            
            fwrite(F_RES1,[pr_coun1;err1;diag(P1).^0.5],'double');
            write_flag1=0;
            
        end
        
        %read the next imu output
        imu_data1=fread(F_IMU1,8,'double');
    end
        
    %INS2 and INS3 computations (lagged state, buffer ins)
    if (pr_coun1>=delay)
        imu2=imu_data2(3:8)-ximu2;
        [Cbn2, vel_n2, pos_n2]=strapdown_ned_dcm_v000(Cbn2, vel_n2, pos_n2, imu2(1:3), imu2(4:6), dt);
        avg_imu2=avg_imu2+imu2;
    
        pr_coun2=pr_coun2+1;
        prp_coun2=prp_coun2+1;
        
        if (mod(prp_coun2,prp_per)==0 && prp_coun2~=0) %propagate the covariance
            %ins 2
            avg_imu2=avg_imu2/prp_coun2;
            [STM2 Q2]=mdl_ned_dcm_v000(pos_n2, vel_n2, Cbn2, avg_imu2(1:3), Aimu, Qimu, Cimu, Rimu, prp_coun2*dt);
            STM2_inv=inv(STM2);
            P2=STM2*P2*STM2'+Q2;
            prp_coun2=0;
            avg_imu2=zeros(6,1);
            
            %update IMU error states
            ximu2=STM2(10:15,10:15)*ximu2;
            
            %ins 2
            dx12_accum=STM2_inv'*dx12_accum;
            dP12_accum=STM2_inv'*dP12_accum*STM2_inv;
            
            %STM
            STM12=STM12*STM2_inv;
            
            write_flag2=1;
        end
        
        while ~isempty(obs_pos2) && pr_coun2==obs_pos2(1)
            %first check wheter the covariance is uptodate
            if (prp_coun2~=0)
                %ins 2
                avg_imu2=avg_imu2/prp_coun2;
                [STM2 Q2]=mdl_ned_dcm_v000(pos_n2, vel_n2, Cbn2, avg_imu2(1:3), Aimu, Qimu, Cimu, Rimu, prp_coun2*dt);
                STM2_inv=inv(STM2);
                P2=STM2*P2*STM2'+Q2;
                prp_coun2=0;
                avg_imu2=zeros(6,1);

                %update IMU error states
                ximu2=STM2(10:15,10:15)*ximu2;

                %ins 2 (FL smoothed) correction values
                dx12_accum=STM2_inv'*dx12_accum;
                dP12_accum=STM2_inv'*dP12_accum*STM2_inv;

                %STM
                STM12=STM12*STM2_inv;
            end
            
            H2=zeros(3,15);
            H2(:,1:3)=eye(3);
            R2=ObsErrDef(1).sR*ObsErrDef(1).sR';
            inno2=pos_n2-obs_pos2(2:4);
                        
            Rinno_inv2=inv(H2*P2*H2'+R2);
            K2=P2*H2'*Rinno_inv2;

            %correct current states
            dx2=K2*inno2;
            pos_n2=pos_n2-dx2(1:3);
            vel_n2=vel_n2-dx2(4:6);
            mx_a=euler2dcm_v000(dx2(7:9));
            Cbn2=mx_a*Cbn2;
            ximu2=ximu2+dx2(10:15);
            P2=(I-K2*H2)*P2*(I-K2*H2)'+K2*R2*K2';

            %update smoother's accumulators
            STM2_inv=inv(I-K2*H2);
            STM12=STM12*STM2_inv;

            dx12_accum=dx12_accum-H2'*Rinno_inv2*inno2;
            dx12_accum=STM2_inv'*dx12_accum;
            
            dP12_accum=dP12_accum-H2'*Rinno_inv2*H2;
            dP12_accum=STM2_inv'*dP12_accum*STM2_inv;
            
            %prepare for the next cycle
            obs_pos2=fread(F_OBS_POS2, 4, 'double');
            write_flag2=1;
            
            if ((pr_coun2==dx_last_up) && (dx_flag==0))     %assign buffer to lagged accumulators to suppress the numerical precision problems (a simple stupid way to reduce the effect of implicit instability of fixed lag smoothers)
                STM12=STM11;
                dx12_accum=dx11_accum;
                dP12_accum=dP11_accum;
                dx_flag=1;
            end

        end
        
       
        %%%write the results
        true_val2 = fread(F_TRU2, 10, 'double');
        imu_err2=fread(F_IMUERR2, 7, 'double');
        if (write_flag2)
            %smoothed values
            dx12=P2*dx12_accum;
            pos_fls=pos_n2-dx12(1:3);
            vel_fls=vel_n2-dx12(4:6);
            mx_a=euler2dcm_v000(dx12(7:9));
            Cbn_fls=mx_a*Cbn2;
            ximu_fls=ximu2+dx12(10:15);
                        
            err2=zeros(15,1);
            err2(1:3)=pos_fls-true_val2(2:4);
            err2(4:6)=vel_fls-true_val2(5:7);
            err2(7:9)=-dcm2euler_v000(Cbn_fls*euler2dcm_v000(true_val2(8:10))');
            err2(10:15)=imu_err2(2:7)-ximu_fls;
            
            %smoothed covariance
            P_fls=(I-P2*dP12_accum)*P2;
            fwrite(F_RES2,[pr_coun2;err2;diag(P_fls).^0.5],'double');
            
            write_flag2=0;
            
%             %write the 2nd ins outputs (not needed at all, just for
%             %checking
%             %read the true values
%             err3=zeros(nstate,1);
%             err3(1:3)=pos2-true_val2(2:4);
%             err3(4:6)=vel2-true_val2(5:7);
%             err3(7:9)=-dcm2euler_v(Cbn2*euler2dcm_v(true_val2(8:10))');
%             err3(10:15)=true_val2(11:16)-imu_err_est2;
%             fwrite(F_RES3,[pr_coun2;err3;diag(P2).^0.5],'double');
        end
        
        %read the next imu output
        imu_data2=fread(F_IMU2,8,'double');
    end
end


fclose all;
figure;
res1=readbin_v000([PATH '\res1.bin'],31);
res2=readbin_v000([PATH '\res2.bin'],31);
plot(res1(1,:)/100,abs(res1(2,:))*6e6)
hold on
plot(res2(1,:)/100,abs(res2(2,:))*6e6,'r')
plot(res1(1,:)/100,(res1(17,:))*6e6,'--');
plot(res2(1,:)/100,(res2(17,:))*6e6,'r--');
title('lat. error (Note that after each update smt=filt for 4 seconds as fls period=11)');
ylabel('m');
xlabel('sec');
legend('filter err.', 'FLS err.','filter std', 'FLS std');
grid;
return;


figure;
yy=6;
plot(res1(1,:)/100,abs(res1(yy,:)))
hold on
plot(res2(1,:)/100,abs(res2(yy,:)),'r')
plot(res1(1,:)/100,(res1(16+yy,:)),'--');
plot(res2(1,:)/100,(res2(16+yy,:)),'r--');
xlabel('sec');
legend('filter', 'fls');
grid;
