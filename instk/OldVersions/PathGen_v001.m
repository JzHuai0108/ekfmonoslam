%%
function rstat=PathGen(dir_name, ini_pva, mot_def, out_typ, sim_mode, sen_conf, imu_conf, vib_param, rstat)

if (rstat(1)~=0)
    randn('state',rstat);
else
    rstat=randn('state');
end

DEBUG=0;    %runs an additional conventional ins to see how compatible the trajectory generator and ins outputs are. (usually required to see the effect numerical approximations under vibration)

%% Frequencies
outfreq=out_typ(1,2);    %imu output frequency
intmul=out_typ(1,1);      %simulation freq multiplier. simulation freq=outfreq*intmul
dt=1/outfreq/intmul;

%%Additional Sensor configurations
nsen=0;     %number of additional sensor units
if (~isempty(sen_conf))    %there exist additional sensors
    nsen=size(sen_conf,1);
    sen_la=zeros(3,nsen);
    sen_or=zeros(3,3,nsen); %dcm from master to sensor
    sen_typ=sen_conf(:,7);
    for in=1:nsen
        sen_la(:,in)=sen_conf(in,1:3); %lever arm in the master body frame
        sen_or(:,:,in)=euler2dcm_v000(sen_conf(in,4:6)); %transformation matrix from sensor to master (euler angles defined in master);
    end
end

%%Additional imu configurations;
nimu=0;     %number of additional imu units
if (~isempty(imu_conf))    %there exist additional sensors
    nimu=size(imu_conf,1);
    imu_la=zeros(3,nimu);
    imu_or=zeros(3,3,nimu); %dcm from master to sensor
    imu_vib=zeros(3,3,nimu);
    for in=1:nimu
        imu_la(:,in)=imu_conf(in,1:3); %lever arm in the master body frame
        imu_or(:,:,in)=euler2dcm_v000(imu_conf(in,4:6)); %transformation matrix from sensor to master (euler angles defined in master);
        imu_vib(:,:,in)=eye(3);     %default value=no vibration
    end
end

%%vibration parameters for additional imus
if (nimu>0)
    vib_state=zeros(3,nimu);
    vib_model=zeros(2,nimu);
    if (~isempty(vib_param))
        for in=1:nimu
            vib_model(:,in)=vib_param(in,1:2)';
            vib_state(3,in)=randn(1)*vib_param(in,3);     %initial vibration angles
            imu_vib(:,:,in)=euler2dcm_v000(vib_state(:,in));
        end
    end
end

%% Path Gen Command Filter and feedback parameters (somehow arbitrary)
INP_FILT_A=eye(3)*0.9;
INP_FILT_B=eye(3)-INP_FILT_A;
att_n_dot=zeros(3,1);
vel_b_dot=zeros(3,1);

%% Open output files
fmimu=fopen([dir_name '\mimu.bin'],'wb');  %sensor outputs
fmnav=fopen([dir_name '\mnav.bin'],'wb');  %master navigation solution

if (nimu>0) %for each additional imus, create additional navigation and imu files
    fsnav=zeros(nimu,1);    %navigation solution files (for each full imu)
    fsimu=zeros(nimu,1);    %imu files
    for in=1:nimu
        fsnav(in)=fopen([dir_name '\snav' int2str(in) '.bin'],'wb');
        fsimu(in)=fopen([dir_name '\simu' int2str(in) '.bin'],'wb');
    end
end

nout=size(out_typ,1);
for in=2:nout    %measurement files
    switch out_typ(in,1)
        case {1}
            fgps=fopen([dir_name '\gps.bin'],'wb');
        case {2}
            fodo=fopen([dir_name '\odo.bin'],'wb');
    end
end

%% convert times into indeces
for in=2:nout;
    out_typ(in,2)=intmul*round(outfreq/out_typ(in,2));
end
for in=1:size(mot_def,1)
    mot_def(in,6)=round(mot_def(in,6)*outfreq*intmul);
end

%%%%%%%%%%% Start Computations %%%%%%%%%%%%%%%
sim_cnt=0;
fflag=1;
acc_tot=zeros(3,1);
gyro_tot=zeros(3,1);
odo_dist=0;

%master navigator states
att_n=ini_pva(:,3);
Cbn=euler2dcm_v000(att_n);
pos_n=ini_pva(:,1);
pos_dn=zeros(3,1);  %total position change
vel_b=ini_pva(:,2);
vel_n=Cbn*vel_b;
[Rn, Re, g, sL, cL, WIE]=geoparam_v000(pos_n);

%Write initial values to the file
fwrite(fmimu,[0;zeros(3,1);zeros(3,1);zeros(nsen,1)],'double');
if (sim_mode==0)
    vr_a=[0;pos_n+pos_dn;vel_b;att_n];
else
    vr_a=[0;pos_dn;vel_b;att_n];
end
fwrite(fmnav, vr_a,'double');

if (nimu>0)
    if (sim_mode==0)
        [spos_n, svel_b, sCsn, simu]=la_imu_v000(zeros(3,1), zeros(3,1), zeros(3,1), pos_n, vel_b, Cbn, imu_la, imu_or, 0, zeros(3,1), imu_vib);
    elseif (sim_mode==1)
        [spos_n, svel_b, sCsn, simu]=la_imu_v000(zeros(3,1), zeros(3,1), zeros(3,1), pos_dn, vel_b, Cbn, imu_la, imu_or, 1, zeros(3,1), imu_vib);
    end
    for in1=1:nimu
        satt_n=dcm2euler_v000(sCsn(:,:,nimu));
        fwrite(fsnav(in1), [0;spos_n(:,in1);svel_b(:,in1); satt_n], 'double');
        fwrite(fsimu(in1), [0;simu], 'double');
    end
end

%%variables for debugging
if (DEBUG)
    %strapdown states (mostly for debugging purposes)
    strap_Cbn=Cbn;
    strap_vel_b=vel_b;
    if (sim_mode==0)
        strap_pos_n=pos_n;
    elseif (sim_mode==1)
        strap_pos_n=pos_dn;
    end

    if (nimu>0)
        strap_spos_n=spos_n;
        for in=1:nimu
            strap_svel_b=svel_b(:,in);
        end
        strap_sCsn=sCsn;
    end
    
    sr_a=sum(mot_def(:,6)); %maximum number of outputs
    deb_vel=zeros(3,sr_a);
    deb_pos=zeros(3,sr_a);
    deb_att=zeros(3,sr_a);
    deb_ctr=0;
end

%Start generation of trajectory
for in=1:size(mot_def,1)
    mmod=mot_def(in,1);
    
    if (mmod==1)
        %Determine the inputs directly
        att_n_dot_com=[mot_def(in,2);mot_def(in,3);mot_def(in,4)];
        vel_b_dot_com=[mot_def(in,5);0;0];
    elseif mmod==2    %new abs. att and vel values to reach
        %determine the commands
        att_n_com=[mot_def(in,2);mot_def(in,3);mot_def(in,4)];
        vel_b_com=[mot_def(in,5);0;0];
    elseif mmod==3    %rel att and vel changes
        att_n_com=[mot_def(in,2)+att_n(1);mot_def(in,3)+att_n(2);mot_def(in,4)+att_n(3)];
        vel_b_com=vel_b+[mot_def(in,5);0;0];
    elseif mmod==4    %abs attitude, rel vel
        att_n_com=[mot_def(in,2);mot_def(in,3);mot_def(in,4)];
        vel_b_com=vel_b+[mot_def(in,5);0;0];
    elseif mmod==5    %rel att, abs vel
        att_n_com=[mot_def(in,2)+att_n(1);mot_def(in,3)+att_n(2);mot_def(in,4)+att_n(3)];
        vel_b_com=[mot_def(in,5);0;0];
    end
        
    att_n_com_filt=att_n;
    vel_b_com_filt=vel_b;

    seg_cnt_lim=sim_cnt+mot_def(in,6);
    seg_next=0;
    
    while ((sim_cnt<seg_cnt_lim) && ~seg_next)
        %determine the inputs
        if (mmod==1)
            %filter the input
            att_n_dot=INP_FILT_A*att_n_dot+INP_FILT_B*att_n_dot_com;
            vel_b_dot=INP_FILT_A*vel_b_dot+INP_FILT_B*vel_b_dot_com;
        else
            %filter the command
            att_n_com_filt=INP_FILT_A*att_n_com_filt+INP_FILT_B*att_n_com;
            vel_b_com_filt=INP_FILT_A*vel_b_com_filt+INP_FILT_B*vel_b_com;

            %close the loop
            att_n_dot=att_n_com_filt-att_n;
            vel_b_dot=vel_b_com_filt-vel_b;
            
            if (norm(att_n_dot)<1e-4 && norm(vel_b_dot)<1e-4)
                seg_next=1;
            end
        end
        
        %Compute the master imu solution
        [acc, gyro, vel_n_dot, pos_n_dot, Cbn_dot]=comp_master(pos_n+pos_dn, vel_b, att_n, Cbn, vel_b_dot, att_n_dot, sim_mode, g);
            
        %update master results
        acc_tot=acc_tot+acc;
        gyro_tot=gyro_tot+gyro;
        pos_dn=pos_dn+pos_n_dot*dt;
        att_n=att_n+att_n_dot*dt;      
        %Cbn=euler2dcm_v000(att_n);
        Cbn=Cbn+Cbn_dot*dt;
        vel_b=vel_b+vel_b_dot*dt;
        vel_n=Cbn*vel_b;
%         vel_n=vel_n+vel_n_dot*dt;
%         vel_b=Cbn'*vel_n;
        
        %%odometer resuts
        odo_vel=vel_b;
        odo_dist=odo_dist+norm(vel_b)*dt;
        
        %update simulation counter
        sim_cnt=sim_cnt+1;
        
        %write the results
        if (mod(sim_cnt,intmul)==0)
            %%average gyro and acc for the output interval
            gyro_avg=gyro_tot/intmul;
            acc_avg=acc_tot/intmul;
           
            if (fflag)  %only once: set gyro derivative to zero, write initial values to the files
                fflag=0;
                gyro_pre=gyro_avg;
                acc_pre=acc_avg;
            end
            
            %%Compute the sensor outputs
            gyro_der=(gyro_avg-gyro_pre)*outfreq;
            if (nsen>0)
                sensor_out=la_sen_v000(acc_avg, gyro_avg, gyro_der, sen_la, sen_or, sen_typ);
            else
                sensor_out=[];
            end
            fwrite(fmimu,[(sim_cnt/intmul);acc_avg;gyro_avg;sensor_out],'double');
            if (sim_mode==0)
                vr_a=[(sim_cnt/intmul);pos_n+pos_dn;vel_b;att_n];
            else
                vr_a=[(sim_cnt/intmul);pos_dn;vel_b;att_n];
            end
            fwrite(fmnav, vr_a,'double');
            
            %Compute the slave imu outputs
            if (nimu>0)
                %update the vibration
                [gyro_vib, vib_state, imu_vib_new]=update_vib(vib_state, vib_model, imu_vib, imu_or, 1/outfreq);
                
                %compute the slave results
                if (sim_mode==0)
                    [spos_n, svel_b, sCsn, simu]=la_imu_v000(acc_avg, gyro_avg, gyro_der, pos_n+pos_dn, vel_b, Cbn, imu_la, imu_or, 0, gyro_vib, imu_vib);
                elseif (sim_mode==1)
                    [spos_n, svel_b, sCsn, simu]=la_imu_v000(acc_avg, gyro_avg, gyro_der, pos_dn, vel_b, Cbn, imu_la, imu_or, 1, gyro_vib, imu_vib);
                end 
                for in1=1:nimu
                    satt_n=dcm2euler_v000(sCsn(:,:,nimu));
                    fwrite(fsnav(in1), [(sim_cnt/intmul);spos_n(:,in1);svel_b(:,in1);satt_n(:,in1)], 'double');
                    fwrite(fsimu(in1), [(sim_cnt/intmul);simu(:,in1)], 'double');
                end
                imu_vib=imu_vib_new;
            end

            %prepare for the next cycle
            gyro_pre=gyro_avg;
            gyro_tot=zeros(3,1);
            acc_tot=zeros(3,1);
            
            %%debug
            if (DEBUG)
                %%run the sptradown
                if (sim_mode==0)
                    [strap_Cbn, strap_vel_b, strap_pos_n]=strapdown_bd_dcm_v000(strap_Cbn, strap_vel_b, strap_pos_n, acc_avg, gyro_avg, 1/outfreq);
                    for in1=1:nimu
                        [strap_sCsn(:,:,in1), strap_svel_b(:,in1), strap_spos_n(:,in1)]=strapdown_bd_dcm_v000(strap_sCsn(:,:,in1), strap_svel_b(:,in1), strap_spos_n(:,in1), simu(1:3,in1), simu(4:6,in1), 1/outfreq);
                    end
                elseif (sim_mode==1)
                    [strap_Cbn, strap_vel_b, strap_pos_n]=strapdown_pln_dcm_v000(strap_Cbn, strap_vel_b, strap_pos_n, acc_avg, gyro_avg, g, 1/outfreq, 1);
                    for in1=1:nimu
                        [strap_sCsn(:,:,in1), strap_svel_b(:,in1), strap_spos_n(:,in1)]=strapdown_pln_dcm_v000(strap_sCsn(:,:,in1), strap_svel_b(:,in1), strap_spos_n(:,in1), simu(1:3,in1), simu(4:6,in1),g, 1/outfreq,1);
                    end
                end
                
                deb_ctr=deb_ctr+1;
%                 deb_vel(:,deb_ctr)=strap_svel_b-svel_b;
%                 deb_pos(:,deb_ctr)=strap_spos_n-spos_n;
%                 deb_att(:,deb_ctr)=dcm2euler_v000(strap_sCsn(:,:,1)*sCsn');
                deb_vel(:,deb_ctr)=strap_vel_b-vel_b;
                deb_pos(:,deb_ctr)=strap_pos_n-pos_dn-pos_n;
                deb_att(:,deb_ctr)=dcm2euler_v000(strap_Cbn(:,:,1)*Cbn');
            end
            
        end
        
        %%measurement files
        if (nout>1)
            for sr_a=2:nout;
                switch out_typ(sr_a,1)
                    case {1}
                        if (mod(sim_cnt,out_typ(sr_a,2))==0)
                            if (sim_mode==0)
                                vr_a=[(sim_cnt/intmul);vel_n;pos_n+pos_dn];
                            elseif(sim_mode==1)
                                vr_a=[(sim_cnt/intmul);vel_n;pos_dn];
                            end
                            fwrite(fgps,vr_a,'double');
                        end
                    case {2}
                        if (mod(sim_cnt,out_typ(sr_a,2))==0)
                            fwrite(fodo, [(sim_cnt/intmul);odo_vel;odo_dist],'double');
                            odo_pos=0;
                        end
                end
            end
        end
    end
end
if (DEBUG)
    deb_vel(:,deb_ctr)
    deb_pos(:,deb_ctr)
    deb_att(:,deb_ctr)
end
fclose('all');

%vel_typ==1 -> b, vel_typ==2 -> n
function [acc, gyro, vel_n_dot, pos_n_dot, Cbn_dot]=comp_master(pos_n, vel_b, att_n, Cbn, vel_b_dot, att_n_dot, styp, g)
%velocity in n
vel_n=Cbn*vel_b;

%%Calculate Spherical Earth Parameters and gravity
if (styp==0)
    [Rn, Re, g, sL, cL, WIE]=geoparam_v000(pos_n);
    Rn_eff=Rn+pos_n(3);
    Re_eff=Re+pos_n(3);
    gravity=[0;0;-g];
elseif (styp==1)
    gravity=[0;0;-g];
end

%%compute w_nb_n ("b" to "n" in "b")
sh=sin(att_n(3));ch=cos(att_n(3));
w_nb_n=zeros(3,1);
w_nb_n(1)=-sh*att_n_dot(2)+Cbn(1,1)*att_n_dot(1);
w_nb_n(2)=ch*att_n_dot(2)+Cbn(2,1)*att_n_dot(1);
w_nb_n(3)=att_n_dot(3)+Cbn(3,1)*att_n_dot(1);

%%compute w_en_n & w_ie_n
if (styp==0)
    w_en_n=zeros(3,1);
    w_en_n(1)=vel_n(2)/(Re_eff);
    w_en_n(2)=-vel_n(1)/(Rn_eff);
    w_en_n(3)=-vel_n(2)*(sL/cL)/(Re_eff);

    w_ie_n=zeros(3,1);
    w_ie_n(1)=+WIE*cL;
    w_ie_n(2)=0;
    w_ie_n(3)=-WIE*sL;
elseif(styp==1)
    w_en_n=zeros(3,1);
    w_ie_n=zeros(3,1);
end

%attitude derivative
Cbn_dot=Cbn*skew(Cbn'*w_nb_n);

%%velocity derivative
vel_n_dot=Cbn*vel_b_dot+skew(w_nb_n)*vel_n;

%%position derivative
if (styp==0)
    pos_n_dot=zeros(3,1);
    pos_n_dot(1)=vel_n(1)/(Rn_eff);
    pos_n_dot(2)=vel_n(2)/(Re_eff)/cL;
    pos_n_dot(3)=-vel_n(3);
elseif (styp==1)
    pos_n_dot=zeros(3,1);
    pos_n_dot(1)=vel_n(1);
    pos_n_dot(2)=vel_n(2);
    pos_n_dot(3)=vel_n(3);
end

%%gyroscope output
gyro=Cbn'*(w_nb_n+w_en_n+w_ie_n);  

% %%Acceleration output
% acc=Cbn'*(vel_ned_dot+(skew(w_en_n)+2*skew(w_ie_n))*vel_ned+gravity);

%%Acceleration output
vr_a=Cbn'*w_ie_n;
acc=vel_b_dot+skew(vr_a+gyro)*vel_b+Cbn'*gravity;



function [gyro_vib, vib_state_new, imu_vib_new]=update_vib(vib_state, vib_model, imu_vib, imu_or, dt)
nimu=size(vib_model,2);
vib_state_new=zeros(3,nimu);
imu_vib_new=zeros(3,3,nimu);
gyro_vib=zeros(3,nimu);
for in=1:nimu
    %compute vibration induced rotation rate in master frame
    Crig=imu_or(:,:,in);    %rigid Csm
    Cvib=imu_vib(:,:,in);   %vibration Cs's , complete rotation Cs'm=Crig*Cvib
    vib_state_new(:,in)=vib_model(1,in)*vib_state(:,in)+[0;0; vib_model(2,in)*randn(1)];
    gyro_vib(:,in)=Crig*(vib_state_new(:,in)-vib_state(:,in))/dt;
    
    %update the vibration 
    rot=Cvib'*(vib_state_new(:,in)-vib_state(:,in));
    rot_norm=norm(rot);
    sr_a=1-(rot_norm^2/6)+(rot_norm^4/120);
    sr_b=(1/2)-(rot_norm^2/24)+(rot_norm^4/720);
    mx_a=eye(3)+sr_a*skew(rot)+sr_b*skew(rot)*skew(rot);
    Cvib=Cvib*mx_a;
    
    imu_vib_new(:,:,in)=Cvib;
end

