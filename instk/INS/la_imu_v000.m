function [spos_n, svel_b, sCsn, simu]=la_imu(acc, gyro, gyro_der, pos_n, vel_b, Cbn, imu_la, imu_or, typ, gyro_vib, imu_vib)

nimu=size(imu_la,2);
spos_n=zeros(3,nimu);
svel_b=zeros(3,nimu);
sCsn=zeros(3,3,nimu);
simu=zeros(6,nimu);

%vibration information is not compulsory
if (isempty(gyro_vib))
    gyro_vib=zeros(3,nimu);
end

if (isempty(imu_vib))
    imu_vib=zeros(3,3,nimu);
    for in=1:nimu
        imu_vib(:,:,in)=eye(3);
    end
end

if (typ==0)
    [Rn, Re, g, sL, cL, WIE]=geo_param(pos_n);
end
for in=1:nimu
    %Attitude
    Crig=imu_or(:,:,in);
    Cvib=imu_vib(:,:,in);
    Csm=Crig*Cvib;
    Csn=Cbn*Csm;
    sCsn(:,:,in)=Csn;
    
    %position
    la_m=imu_la(:,in);
    if (typ==0)
        spos_n(:,in)=pos_n+[1/(Rn+pos_n(3)) 0 0;0 1/cL/(Re+pos_n(3)) 0; 0 0 -1]*Cbn*la_m;
    elseif (typ==1)
        spos_n(:,in)=pos_n+Cbn*la_m;
    end
    
    %Velocity
    if (typ==0)
        wie_m=Cbn'*[WIE*cL; 0; -WIE*sL];
    else
        wie_m=zeros(3,1);
    end
    svel_b(:,in)=Csm'*(vel_b+cross(gyro-wie_m,la_m));

    %Acceleration
    simu(1:3,in)=Csm'*(acc+cross(gyro,cross(gyro,la_m))+cross(gyro_der,la_m));

    %Rotation rate
    simu(4:6,in)=Csm'*(gyro+gyro_vib(:,in));
end