fclose all;
clear;

TIME.dt=1/400;  %I use actual time difference as dt. an approximate estimate is enough for this.
TIME.npt=4; %Number of subminor intervals
TIME.feedInterval=10; %Time interval between 2 navigation error feedback. (I do not feedback imu errors at all)

%Observation params
OBS.odoR=0.2^2;  %(m/s)^2
OBS.useNonHolo=1;   %0:Don't use, 1:use
OBS.nonHoloR=0.1^2; %(m/s)^2
OBS.gpsR=0.1^2'; %(m)^2
OBS.zuptR=0.1^2;
OBS.LA=[0.2;0;0];

%open files
fimu=fopen('/home/yigiter/Desktop/Dropbox/INS/SAGE/mfiles/data/imudata.bin','r');
fnavres=fopen('/home/yigiter/Desktop/Dropbox/INS/SAGE/mfiles/data/navresult.bin','r');
fcovres=fopen('/home/yigiter/Desktop/Dropbox/INS/SAGE/mfiles/data/covresult.bin','r');

%Odo and gps data are text files. Therefore, it is easier to read them all
%at once. (We do not know how many bytes there are in a row in the text file)
ododata=load('/home/yigiter/Desktop/Dropbox/INS/SAGE/mfiles/data/odo.txt');
gpsdata=load('/home/yigiter/Desktop/Dropbox/INS/SAGE/mfiles/data/gps.txt');

fsmtnavres=fopen('/home/yigiter/Desktop/Dropbox/INS/SAGE/mfiles/data/smtnavresult.bin','w');
fsmtcovres=fopen('/home/yigiter/Desktop/Dropbox/INS/SAGE/mfiles/data/smtcovresult.bin','w');
fbcknavres=fopen('/home/yigiter/Desktop/Dropbox/INS/SAGE/mfiles/data/bcknavresult.bin','w');
fbckcovres=fopen('/home/yigiter/Desktop/Dropbox/INS/SAGE/mfiles/data/bckcovresult.bin','w');
%fdebug=fopen('/home/yigiter/Desktop/Dropbox/INS/SAGE/mfiles/data/debug.bin','w');

%Read initial (nominal) values.
%In fact, you can use whatever initial value you wish. Thereforetically,
%it should not affect the results .In other words, the smoother results
%theoretically invariant to the initial state conditions. However, in
%practice it is not. That is why I use the final filtering value as the
%initial value.
%On the other hand, I used 0 initial condition for imu error states. This
%is because, I am going to propagate the error states as z=Pinv*dx which
%does not require the inverse of STM. If I want to use a non-zero initial
%value for imu error states, I have to explicitly propagate them using the
%inverse differential equations. You should note that navigation error
%states are explicitly propagated in the backwards ins when they are feeded
%back. 

vr_a=fread(fnavres,[30,1],'double');
TIME.interval(1)=vr_a(1);   %Starting time
inillh=[dcm2llh_v000(reshape(vr_a(2:10),[3,3]));vr_a(11)]; %initial position

fseek(fnavres, -30*8,'eof');
fwdnav=fread(fnavres,[30,1],'double');
fwdtime=fwdnav(1);
TIME.interval(2)=fwdnav(1); %stop time. (do not need this)
Cen=reshape(fwdnav(2:10),[3,3]);
height=fwdnav(11); %elipsoid height
Vn=fwdnav(12:14);
qbn=fwdnav(15:18);

%The z=PbvkInv*dx
dz=zeros(21,1);
%Note that this implies dx=0. Therefore, I do not use the imu error
%estimates as the initial nominal errors. It can be done, however, it
%requires that these nominal values must be propagated backwards explicitly
%which is something that I avoid.

%Forward covariance
fseek(fcovres, -232*8, 'eof');
fwdcov=fread(fcovres,[232,1],'double');
Pfwd=mat2vec_v000(fwdcov(2:end));

%inverse of the covaraince
PbckInv=zeros(21);
covupt_time=fwdtime;

% %%debug
% dx=zeros(21,1);
% Pbck=diag([10000 10000 10000 0.1 0.1 0.1 0.001 0.001 0.001 9e-6 9e-6 9e-6 2.3e-11 2.3e-11 2.3e-11 0.01 0.01 0.01 0.01 0.01 0.01]);
% PbckInv=inv(Pbck);


%Adjust the imu file pointer
fseek(fimu,-7*8,'eof');
imudata=fread(fimu,[7,1],'double');
if (imudata(1)>fwdtime)
    fseek(fimu, -7*floor((imudata(1)-fwdtime)/TIME.dt)*8,'cof');
end
imudata=fread(fimu,[7,1],'double');
while(imudata(1)>fwdtime)
    fseek(fimu,-7*8*2,'cof');
    imudata=fread(fimu,[7,1],'double');
end
fseek(fimu,(TIME.npt-1)*7*8,'cof');

%Adjust the Observation counters
odoctr=find(ododata(:,1)>=(fwdtime),1)-1;
gpsctr=find(gpsdata(:,1)>=(fwdtime),1)-1;
if (isempty(gpsctr))
    gpsctr=size(gpsdata,1);
end

%%%Start the backward computations
%Prepare for the loop
imuaccum=zeros(6,1);
curimutime=imudata(1);
lastfeedback=curimutime;
dt=TIME.dt*TIME.npt; %actual value will be computed when the imu is read again.
initime=curimutime;

while (imudata(1)>TIME.interval(1))
    if (initime-imudata(1))>60
        disp(['Process Time:' num2str(imudata(1)/60)]);
        initime=imudata(1);
    end
    %%%%%Propagate the covariance and states
    %ododata(odoctr,1)=0;
    if (curimutime<covupt_time) && ((curimutime<=fwdtime) || curimutime<(ododata(odoctr,1)+dt) || curimutime<(gpsdata(gpsctr,1)+dt))
        %%update the covariance
        covdt=covupt_time-curimutime;
        [STM, Qd]=sys_metric_phipsi_v000(Cen, height, Vn, qbn, imuaccum(4:6)/covdt, imuaccum(1:3)/covdt, covdt,2,1); %Note:this is a forwards model as it must be.
        mx_a=STM'*inv(eye(21)+PbckInv*Qd);
        PbckInv=mx_a*PbckInv*STM;
        
        %%Update the (inverse covariance weighted) state estimates
        dz=mx_a*dz;

%         %%Debug (assuming that we have a Pbck)
%         mx_a=inv(STM);
%         Pbck=mx_a*Pbck*mx_a'+mx_a*Qd*mx_a';
%         dx=mx_a*dx;
%         mx_c=inv(PbckInv);
%         mx_b=Pbck-mx_c;
%         vr_a=dx-mx_c*dz;
        
        imuaccum=zeros(6,1);
        covupt_time=curimutime;
        
        %Apply ZUPT.
        if ((curimutime>12125 && curimutime<12174) || (curimutime>17250 && curimutime<17337))
            inno=Vn;
            H=[zeros(3) eye(3) zeros(3,15)];
            R=eye(3)*OBS.zuptR;
            
            %Kalman
            mx_a=H'*inv(R);
            PbckInv=PbckInv+mx_a*H;
            dz=dz+mx_a*inno;
            
%         %%Debug
%         K=Pbck*H'*inv(H*Pbck*H'+R);
%         dx=dx+K*(inno-H*dx);
%         Pbck=(eye(21)-K*H)*Pbck*(eye(21)-K*H)'+K*R*K';
%         mx_c=inv(PbckInv);
%         mx_b=Pbck-mx_c;
%         vr_a=dx-mx_c*dz;
        end
    end
    
    %%%%Apply the observations
    %Note: Although I do not feedback errors (not all of them at least),
    %I do not use dz during the innovation computations.This is because update equation
    %of the information form also takes care of this automatically.
    
    %%odometer updates
    if (curimutime<(ododata(odoctr,1)+dt))  %we have an odometer data
        Cnb=quat2dcm_v000(qbn)';
        if (OBS.useNonHolo)
            inno=[ododata(odoctr,2);0;0]-Cnb*Vn;
            H=[zeros(3) -Cnb Cnb*skew(Vn) zeros(3,12)];
            R=diag([OBS.odoR, OBS.nonHoloR, OBS.nonHoloR]);
        else
            inno=ododata(odoctr,2)-Cnb(1,:)*Vn;
            H=[zeros(1,3) -Cnb(1,:) Cnb(1,:)*skew(Vn) zeros(1,12)];
            R=OBS.odoR;
        end
        
        %Kalman
        mx_a=H'*inv(R);
        PbckInv=PbckInv+mx_a*H;
        dz=dz+mx_a*inno;

%         %%Debug
%         K=Pbck*H'*inv(H*Pbck*H'+R);
%         dx=dx+K*(inno-H*dx);
%         Pbck=(eye(21)-K*H)*Pbck*(eye(21)-K*H)'+K*R*K';
%         mx_c=inv(PbckInv);
%         mx_b=Pbck-mx_c;
%         vr_a=dx-mx_c*dz;
        
        %new odo data index
        odoctr=odoctr-1;
        
        if (odoctr==0)
            odoctr=1;
            ododata(odoctr,1)=-1;
        end
    end
    
    %%gps updates
    if ((curimutime<=gpsdata(gpsctr,1)+dt))
%         Rappx=6378137/(sqrt(1.0-0.00669437999014*Cen(3,3)^2))+height;
%         gpssL=sin(gpsdata(gpsctr,2)*pi/180);
%         gpscL=cos(gpsdata(gpsctr,2)*pi/180);
%         gpssl=sin(gpsdata(gpsctr,3)*pi/180);
%         gpscl=cos(gpsdata(gpsctr,3)*pi/180);
%         vr_a=Rappx*Cen(1:2,:)*[gpscL*gpscl;gpscL*gpssl;gpssL];
%         vr_b=quatrot_v000(qbn,OBS.LA,0);
%         inno=[-vr_a;gpsdata(gpsctr,4)-height]+vr_b;
        
        gpsecef=ecef2geo_v000([gpsdata(gpsctr,2)/180*pi;gpsdata(gpsctr,3)/180*pi;gpsdata(gpsctr,4)],1);
        Re=6378137/(sqrt(1.0-0.00669437999014*Cen(3,3)^2));
        posecef=-[(Re+height)*Cen(3,1);(Re+height)*Cen(3,2);(Re*(1-0.00669437999014)+height)*Cen(3,3)];
        inno=Cen*(posecef-gpsecef)+quatrot_v000(qbn,OBS.LA,0);

        disp(inno);
        
        H=[eye(3) zeros(3,18)];
        R=eye(3)*OBS.gpsR;
        
        %Kalman
        mx_a=H'*inv(R);
        PbckInv=PbckInv+mx_a*H;
        dz=dz+mx_a*inno;
        
%         %%Debug
%         K=Pbck*H'*inv(H*Pbck*H'+R);
%         dx=dx+K*(inno-H*dx);
%         Pbck=(eye(21)-K*H)*Pbck*(eye(21)-K*H)'+K*R*K';
%         mx_c=inv(PbckInv);
%         mx_b=Pbck-mx_c;
%         vr_a=dx-mx_c*dz;

        gpsctr=gpsctr-1;
        if (gpsctr==0)
            gpsctr=1;
            gpsdata(gpsctr,1)=-1;
        end
    end
    
    
    %Combine and write solutions
    if (curimutime<=fwdtime)
        %Write backward solution
        fwrite(fbcknavres,[imudata(1);Cen(:);height;Vn;qbn;dz],'double');
        fwrite(fbckcovres,[imudata(1);mat2vec_v000(PbckInv)],'double');
        
%         %%Debug
%         [qbn_bck, Vn_bck, Cen_bck, h_bck]=correctnav_Cen_v000(qbn, Vn, Cen, height, dx(1:9), 2, 2);
%         fwrite(fdebug,[imudata(1);Cen_bck(:);h_bck;Vn_bck;qbn_bck;...
%         Cen(:);height;Vn;qbn],'double');
        
        %Compute the smoothed solution
        [P_smt, Cen_smt, h_smt, Vn_smt, qbn_smt, imu_smt, dxDif]=combine_2filt_v001(fwdnav(2:end), fwdcov(2:end), Cen, height, Vn, qbn, dz, PbckInv);       
        
        %Write the smoothed solution
        fwrite(fsmtnavres,[imudata(1);Cen_smt(:);h_smt;Vn_smt;qbn_smt;imu_smt],'double');
        fwrite(fsmtcovres,[imudata(1);mat2vec_v000(P_smt)],'double');

        %Read the next forward solution
        fseek(fnavres, 2*-30*8,'cof');
        fwdnav=fread(fnavres,[30,1],'double');
        fwdtime=fwdnav(1);
        
        fseek(fcovres, 2*-232*8, 'cof');
        fwdcov=fread(fcovres,[232,1],'double');
    end
    
    %Feedback the navigation errors
    if ((lastfeedback-curimutime)>TIME.feedInterval)
        %We cannot provide feedback based on dz because there is no
        %guarantee that PbckInv will ever be invertible.
        %One way might be to analyze the svd's of PbckInv and compute some
        %approximate inverse values for navigation states. However, due to
        %the magnitude difference between the nav and imu states, this method
        %is also quite problematic.
        
        %therefore, here I use a simple shortcut. I use smoother error
        %estimates as the feedback values. (We can thereoretically
        %add/subctract any value from estimates and nominal values at the
        %same time without affecting the optimality)
        
        %Use smoothing error values to correct the nav states
        [qbn, Vn, Cen, height]=correctnav_Cen_v000(qbn, Vn, Cen, height, dxDif(1:9), 2, 2);
        %Update dz accordingly
        dz=dz-PbckInv(:,1:9)*dxDif(1:9);
        
        %Note that I do not change the imu error states above. (This is
        %because I intentionally avoid implementing a nominal backwards imu
        %model)
        
        lastfeedback=curimutime;
    end
    
    %Read the next imu
    preimutime=curimutime;
    fseek(fimu, -2*TIME.npt*7*8,'cof');
    imudata=fread(fimu, [7,TIME.npt], 'double');
    if (isempty(imudata))
        break;
    end
    curimutime=imudata(1,1);
    dt=preimutime-curimutime;
    
    %%%process the imuda data backwards
    %Simply accumulate the 4 data samples
    angleinc=sum(imudata(2:4,:),2);
    velinc=sum(imudata(5:7,:),2);
    
    %Note: there is no imu compansation as in the forward filter. Because,
    % I am not going to feedback the imu errors.
    imuaccum=imuaccum+[angleinc;velinc];
    
    %run strapdown backwards
    [qbn, Vn, Cen, height]=strapdown_bckwrd_Cen_quat_v000(qbn, Vn, Cen, height, velinc, angleinc, dt);
    
end

fclose all;


%%%%%%%%%%%%%%%SHOW RESULTS
fff=readbin_v000('/home/yigiter/Desktop/Dropbox/INS/SAGE/mfiles/data/navresult.bin',30);
bbb=readbin_v000('/home/yigiter/Desktop/Dropbox/INS/SAGE/mfiles/data/bcknavresult.bin',18+21);
sss=readbin_v000('/home/yigiter/Desktop/Dropbox/INS/SAGE/mfiles/data/smtnavresult.bin',30);
zs=posdiff_v000([sss([4,7,10],:)]',sss(11,:)', inillh);
zf=posdiff_v000([fff([4,7,10],:)]',fff(11,:)', inillh);
zb=posdiff_v000([bbb([4,7,10],:)]',bbb(11,:)', inillh);
%zd=posdiff_v000([ddd([4,7,10],:)]',ddd(11,:)', inillh);
zg=posdiff_v000(gpsdata(:,2:3)*pi/180,gpsdata(:,4), inillh);

plot(zf(:,1), zf(:,2))
hold on
plot(zg(:,1), zg(:,2),'r+')
plot(zs(:,1), zs(:,2),'r')
plot(zb(:,1), zb(:,2),'k')
%plot(zd(:,1), zd(:,2),'--r')

% zs1=zs;
% zf1=zf;
% zb1=zb;
% sss1=sss;
% fff1=fff;
% save('benisil1.mat','zs1','zf1','zb1','sss1','fff1');
