%nrow: 1 + number of rows in the IMU file (For an 6dof IMU this must be 7)
%act_row: specify which sensors (i.e. rows) to add error.
%errdefs: error models for the sensors specified by the act_row
%tem_mod: temperature model (characterize the change of temperature, not the sensor's temp characteristics)
%rstat: seed for the number. if it is empty, it have the generator shuffled.

function [rstat]=AddIMUErr(finp_name, fout_name, ferr_name, nrow, act_row, errdefs, tem_mod, rstat)

%randomize
% if (isempty(rstat))
%     rstat=get(RandStream.getDefaultStream,'State');
% else
%     set(RandStream.getDefaultStream,'State',rstat);
% end
if (isempty(rstat))
    rstat=rng('shuffle');
else
    rng(rstat);
end

if isempty(tem_mod) %if no tem_mod is specified, use constant temperature=0.
    tem_mod.A=1;
    tem_mod.B=0;
    tem_mod.u=0;
end

%input data
if (~isempty(finp_name))
    fimu=fopen(finp_name,'rb');
    imu_inp=fread(fimu,[nrow inf],'double');
    fclose(fimu);
    imu_inp=imu_inp(act_row,:);
else    %if there is no input file, use constant values as input (act_row=inp values, nrow=#inp)
    imu_inp=act_row*ones(1,nrow);
end

ndat=size(imu_inp,2);

%output files
feimu=fopen(fout_name,'wb');
ferr=fopen(ferr_name,'wb');

if (isempty(errdefs))   %do not add error
     out=[0:ndat-1;zeros(1,ndat);imu_inp];
     fwrite(feimu,out,'double');
else    %generate and add errors
    %Generate the system model for the TI imu errors
    [Ati, Bti, Cti, Dti, sPti]=imu_modTI_v000(errdefs);
    nst_ti=size(Ati,1);
    st_ti=sPti*randn(nst_ti,1);

    tem=0;
    tem_dif=0;
    [Atv, Btv, Ctv, sPtv]=imu_modTV_v000(errdefs, tem_dif, tem, 1);
    if (~isempty(Atv))
        nst_tv=size(Atv,1);
        st_tv=sPtv*randn(nst_tv,1);
    else
        nst_tv=0;
    end

    %add errors
    nsen=size(Cti,1);
    if (nst_tv==0) %only ti part
        for in=1:ndat
            imu_out=imu_inp(:,in)+Cti*st_ti+Dti*randn(nsen,1);
            fwrite(ferr,[in-1;st_ti],'double');
            fwrite(feimu,[in-1;tem;imu_out],'double');

            %new error values
            tem=tem_mod.A*tem+tem_mod.B*randn(1)+tem_mod.u;
            st_ti=Ati*st_ti+Bti*randn(nst_ti,1);
        end
    else    %both ti and tv parts are generated
        for in=1:ndat
            imu_out=imu_inp(:,in)+Cti*st_ti+Dti*randn(nsen,1)+Ctv*st_tv;
            fwrite(ferr,[in-1;st_ti;st_tv],'double');
            fwrite(feimu,[in-1;tem;imu_out],'double');

            %new error values
            sr_a=tem_mod.A*tem+tem_mod.B*randn(1)+tem_mod.u;
            tem_dif=sr_a-tem;
            tem=sr_a;
            [Atv, Btv, Ctv]=imu_modTV_v000(errdefs, tem_dif, tem, 0);
            st_tv=Atv*st_tv+Btv*randn(nst_tv,1);
            st_ti=Ati*st_ti+Bti*randn(nst_ti,1);
        end
    end
end
fclose(feimu);
fclose(ferr);
