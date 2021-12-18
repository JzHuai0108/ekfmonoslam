%Standard methods.
%Cbn_est=(I-S(e))Cbn, e=E[acc_err;ref_err]
%%TODO: The argument names became meaningless. I need to revise this.

function [Cnb E]=alingc_std(acc, ref, Llh, mtd)
    E=[];
    if (mtd==0) %gyroscompass
        %normalize acc and gyro;
        gyro=ref;
        acc=acc/norm(acc);
        gyro=gyro/norm(gyro);
        cL=cos(Llh(1));
        sL=cos(Llh(2));
        Cnb=[cross(acc,gyro) gyro acc]*[0 -1/cL 0;1/cL 0 0;-sL/cL 0 -1];
    elseif (mtd==1) %Rogers, Titterton etc (gyrocompass)
        [Rn, Re, g, sL, cL, WIE_E]=geoparam_v000(Llh);
        gyro=ref;
        Cnb=[cross(acc,gyro) gyro acc]*[0 -1/cL/WIE_E/g 0;1/cL/WIE_E 0 0;-sL/cL/g 0 -1/g];
    elseif (mtd==2) %ref is any vector in the x-z plane (could be gyro or compass output (ignoring declination or aligning to mag north))
        nref=norm(ref);
        nacc=norm(acc);
        vr_z=-acc/nacc;
        vr_y=cross(vr_z,ref/nref);
        nvr_y=norm(vr_y);
        vr_y=vr_y/nvr_y;
        vr_x=cross(vr_y,-acc/nacc);

        Cnb=[vr_x vr_y vr_z];

        %Error matrix (completely ignores the effect of norm errors. (hence, not
        %so correct)
        E=zeros(3,6);
        E(1,1:3)=(Cnb(:,2)/nacc)';
        E(2,1:3)=(-Cnb(:,1)/nacc)';
        E(3,4:6)=(cross(Cnb(:,1),acc/nacc)/nref/nvr_y)';
        E(3,1:3)=(cross(Cnb(:,1),-ref/nref)/nacc/nvr_y)';
    end
end
