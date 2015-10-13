
% the Cb2imu can be determined by considering that the body frame is
% related to the n-frame by only a rotation on the down axis, and by coarse
% calibration, we can determine Cmem2n, then Cmem2b and Cn2mem
% this function requires its data are collected in static mode
function qs2b=getqimu2body(imufile, inillh, startTime, dt)
alignEpoch=40;
preimudata=CQueue();    % record the previous imu data
fimu=fopen(imufile,'r');
fgetl(fimu);            % remove the header
hstream= fgetl(fimu);
mass=textscan(hstream,'%f','delimiter',',');
imudata=mass{1};
imudata(2:7,1)=imudata(2:7,1)*dt;
while(imudata(7,1)==0||imudata(1,1)<startTime)
    preimudata.push(imudata(:,end));%record previous imudata
    if(preimudata.size()>alignEpoch)
        preimudata.pop();
    end
    hstream= fgetl(fimu);
    mass=textscan(hstream,'%f','delimiter',',');
    imudata=mass{1};
    imudata(2:7,1)=imudata(2:7,1)*dt;
end
fclose(fimu);
databuf=reshape(cell2mat(preimudata.content()),7,preimudata.size());
%the method of Yuksel does not work out as well as Yudan's coarse
%alignment
%     meandat=mean(databuf(2:end,:),2)/options.dt;
%     %Compute the alignment between the geodetic and the body frame
%     [Fc wen_n wie_n g]=geoparam_v001(2, Cen(:,3), height, zeros(3,1));
%     [Cnb E err]=align_opt_v000([meandat(4:6), meandat(1:3)],[[0;0;-g], wie_n] ,[1 0]); %Note:E defined for errors on imu. Therefore, the order of vectors is important
%     qbn=cbn2quat(Cnb');
% nav(8:10) corresponds to attitude, roll, pitch, yaw of Cs2n
navdata = coarse_alignment(databuf',inillh, [],0, 3);
qs2b=att2qua([navdata(8:9),0]);% split the rotation vector
end