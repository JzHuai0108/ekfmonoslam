function M=cp_M(imu_conf)
nsen=size(imu_conf,1);
M=zeros(nsen,6);

%Form M assuming full col rank.
for in=1:nsen
    vr_a=imu_conf(in,4:6);
    vr_a=vr_a/norm(vr_a);
    if (imu_conf(in,7)==1) %accelerometer
        M(in,1:3)=vr_a;
    else %gyroscope
        M(in,4:6)=vr_a;
    end
end

%remove all-zero columns
nzr_col=[];
for in=1:6
    if (sum(M(:,in)==0)==nsen)
        %disp('cp_M: Zero col');
    else
        nzr_col=[nzr_col,in];
    end
end

M=M(:,nzr_col);