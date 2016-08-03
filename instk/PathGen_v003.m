%Simplified version of PathGen_v002 for only a single INS and without any vibration model.

%dir_name : Directory to write outputs (must include trailing slash)
%ini_pva :[pos_in_nav(lLh), vel_in_body, body_nav_euler];
%out_typ :[[simulation_freq_multiplier imu_output_freq];[1 gps_freq]; [2 odo_freq]]
%sim_mode :
%   0:nav frame simulation,
%   1:simulation in a planar frame with constant g (position output in meters) (ini_pos must also be in meters)
%mot_def : Definition of motions. Each row defines a segment.
%   mot_def(:,1)=motion type (see the code for the options)
%   mot_def(:,2:5)=parameters of the motion (see the code)
%   mot_def(:,6)=maximum allowed time for the given segment

%Output Files:
%mimu.bin: imu_data=[out_ind, acc, rotation_rate] (time=out_ind/imu_freq)
%mnav.bin: true_nav_data=[out_ind, pos_n or pos_dn, vel_b;att_n] (sim_mode determines position mechanization)

function PathGen(dir_name, ini_pva, mot_def, out_typ, sim_mode)

%Run an additional INS to see how compatible the trajectory generator and ins outputs are.
%(This is usually required to see the effect of numerical approximations under vibration)
DEBUG=0; %0 --> Don't use an extra INS, 1 --> Use additional INS for debugging

% Frequencies
outfreq=out_typ(1,2);       %imu output frequency
intmul=out_typ(1,1);        %simulation freq multiplier. simulation freq=outfreq*intmul
simfreq=outfreq*intmul;     %simulation frequency
dt=1/simfreq;

% Path Gen Command Filter and feedback parameters (somehow arbitrary. Just made them up)
time_const=1; %Response time constant (in sec)
INP_FILT_A=eye(3)*exp(-1/time_const/simfreq);
INP_FILT_B=eye(3)-INP_FILT_A;
att_n_dot=zeros(3,1);
vel_b_dot=zeros(3,1);

% Open output files
fmimu=fopen([dir_name 'mimu.bin'],'wb');  %sensor outputs
fmnav=fopen([dir_name 'mnav.bin'],'wb');  %master navigation solution

% Measurement files
nout=size(out_typ,1);
for in=2:nout
    switch out_typ(in,1)
        case {1}
            fgps=fopen([dir_name 'gps.bin'],'wb');
        case {2}
            fodo=fopen([dir_name 'odo.bin'],'wb');
    end
end

% Convert times into indices
for in=2:nout;
    out_typ(in,2)=round(simfreq/out_typ(in,2));
end
for in=1:size(mot_def,1)
    mot_def(in,6)=round(mot_def(in,6)*simfreq);
end

%%%%%%%%%%% Start Computations %%%%%%%%%%%%%%%
sim_cnt=0;
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
fwrite(fmimu,[0;zeros(3,1);zeros(3,1)], 'double');
if (sim_mode==0)
    vr_a=[0;pos_n+pos_dn;vel_b;att_n];
else
    vr_a=[0;pos_dn;vel_b;att_n];
end
fwrite(fmnav, vr_a,'double');

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
    elseif mmod==6    %bank to turn (before using this, be sure that roll and pitch is zero)
        % bank mot_def(in,2) radians to reach mot_def(in,4) radians azimuth change.
        % if mot_def(in,4)==0, only roll angle is used to determine the exit conditon
        %                      (i.e. use [6 0 0 0 0 max_time] to return to level flight from a banked manevour)
        % abs roll, rel heading
        % segment ends when required azimuth is reached
        att_n_com=[mot_def(in,2); 0; mot_def(in,4)+att_n(3)];
        vel_b_com=[0;0;0];
    end

    seg_cnt_lim=sim_cnt+mot_def(in,6);
    seg_next=0;
    
    while ((sim_cnt<seg_cnt_lim) && ~seg_next)
        %determine the inputs
        if (mmod==1) %control is based on att and vel derivative
            %filter the input
            att_n_dot=INP_FILT_A*att_n_dot+INP_FILT_B*att_n_dot_com;
            vel_b_dot=INP_FILT_A*vel_b_dot+INP_FILT_B*vel_b_dot_com;
        elseif mmod>1 && mmod<6 %control is based on att and vel values
            %Generate the command
            att_err=att_n_com-att_n;
            vel_err=vel_b_com-vel_b;
            
            %filter the command
            att_n_dot=INP_FILT_A*att_n_dot+INP_FILT_B*att_err;
            vel_b_dot=INP_FILT_A*vel_b_dot+INP_FILT_B*vel_err;
            
            if (norm(att_err)<1e-2 && norm(vel_err)<1e-2 && norm(att_n_dot)<1e-4 && norm(vel_n_dot)<1e-4)   
                seg_next=1; %when we reach the segment objective, continue to the next segment.
            end
        elseif mmod==6 %bank to turn for aircrafts
            %For this to be meaningful, pitch must be zero.
            roll_err=att_n_com(1)-att_n(1);
            
            %filter the command
            att_n_dot(1)=INP_FILT_A(1,1)*att_n_dot(1)+INP_FILT_B(1,1)*roll_err;
            att_n_dot(2)=0;
            att_n_dot(3)=g*tan(att_n(1))/norm(vel_b);
            
            vel_b_dot=zeros(3,1);
            
            if (mot_def(in,4)==0)   %is there any heading change specification?
                att_err=att_n_com(1)-att_n(1);  %NO: Use only roll angle to control (for levelling)
                if (norm(att_err)<1e-2 && norm(att_n_dot)<1e-4)   
                    seg_next=1; %when we reach the segment objective, continue to the next segment.
                end
            else
                att_err=att_n_com(3)-att_n(3);  %Yes: Use only heading angle to control (for bank to turn)
                
                if (norm(att_err)<1e-2)   
                    seg_next=1; %when we reach the segment objective, continue to the next segment.
                end
            end
        end
        
        %Compute the master imu solution
        [acc, gyro, vel_n_dot, pos_n_dot, Cbn_dot]=comp_master(pos_n+pos_dn, vel_b, att_n, Cbn, vel_b_dot, att_n_dot, sim_mode, g);
            
        %update master results
        acc_tot=acc_tot+acc;
        gyro_tot=gyro_tot+gyro;
        pos_dn=pos_dn+pos_n_dot*dt;
        att_n=att_n+att_n_dot*dt;      
        Cbn=euler2dcm_v000(att_n);
        %Cbn=Cbn+Cbn_dot*dt;
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
            %average gyro and acc for the output interval
            gyro_avg=gyro_tot/intmul;
            acc_avg=acc_tot/intmul;
           
            %IMU outputs
            fwrite(fmimu,[(sim_cnt/intmul);acc_avg;gyro_avg],'double');
            
            %NAV outputs
            if (sim_mode==0)
                vr_a=[(sim_cnt/intmul);pos_n+pos_dn;vel_b;att_n];
            elseif (sim_mode==1)
                vr_a=[(sim_cnt/intmul);pos_dn;vel_b;att_n];
            end
            fwrite(fmnav, vr_a,'double');

            %prepare for the next cycle
            gyro_tot=zeros(3,1);
            acc_tot=zeros(3,1);
            
            %%debug
            if (DEBUG)
                %%run the sptradown
                if (sim_mode==0)
                    [strap_Cbn, strap_vel_b, strap_pos_n]=strapdown_bd_dcm_v000(strap_Cbn, strap_vel_b, strap_pos_n, acc_avg, gyro_avg, 1/outfreq);
                elseif (sim_mode==1)
                    [strap_Cbn, strap_vel_b, strap_pos_n]=strapdown_pln_dcm_v000(strap_Cbn, strap_vel_b, strap_pos_n, acc_avg, gyro_avg, g, 1/outfreq, 1);
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
                    case {1} %GPS
                        if (mod(sim_cnt,out_typ(sr_a,2))==0)
                            if (sim_mode==0)
                                vr_a=[(sim_cnt/intmul);vel_n;pos_n+pos_dn];
                            elseif(sim_mode==1)
                                vr_a=[(sim_cnt/intmul);vel_n;pos_dn];
                            end
                            fwrite(fgps,vr_a,'double');
                        end
                    case {2} %ODO
                        if (mod(sim_cnt,out_typ(sr_a,2))==0)
                            fwrite(fodo, [(sim_cnt/intmul);odo_vel;odo_dist],'double');
                        end
                end
            end
        end
    end
    
    if (DEBUG)
        disp(['Segment ' num2str(in) ' ends: ' num2str(sim_cnt/simfreq) 'Sec']);
    end
end
if (DEBUG)
    deb_vel(:,deb_ctr)
    deb_pos(:,deb_ctr)
    deb_att(:,deb_ctr)
end
fclose('all');


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

%compute w_nb_n (rotation rate of "b" wrt "n" defined in "n")
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
elseif(styp==1) %navigation frame is fixed
    w_en_n=zeros(3,1);
    w_ie_n=zeros(3,1);
end

%attitude derivative
Cbn_dot=Cbn*skew(Cbn'*w_nb_n);

%%velocity derivative
vel_n_dot=Cbn*vel_b_dot+skew(w_nb_n)*vel_n;

%%position derivative
if (styp==0) %lat, long, height derivative)
    pos_n_dot=zeros(3,1);
    pos_n_dot(1)=vel_n(1)/(Rn_eff);
    pos_n_dot(2)=vel_n(2)/(Re_eff)/cL;
    pos_n_dot(3)=-vel_n(3);
elseif (styp==1) %position in a cartesian 'n' frame 
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
