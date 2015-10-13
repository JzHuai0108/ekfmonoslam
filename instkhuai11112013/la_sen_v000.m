function sensor_out=la_sen(acc, gyro, gyro_der, sen_la, sen_or, sen_typ)
nsen=size(sen_typ,1);
sensor_out=zeros(nsen,1);
for in=1:nsen;
    %lever arm
    Csm=sen_or(:,:,in);
    la=sen_la(:,in);
    
    %Acceleration
    acc_sen=Csm'*(acc+cross(gyro,cross(gyro,la))+cross(gyro_der,la));

    %Rotation rate
    gyro_sen=Csm'*gyro;

    if (sen_typ(in)==1)
        sensor_out(in)=acc_sen(1);
    elseif (sen_typ(in)==2)
        sensor_out(in)=gyro_sen(1);
    end
end