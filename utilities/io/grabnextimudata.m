function [fimu, imudata]=grabnextimudata(fimu, preimutime, imuFileType)
% grab the next imu data that has a timestamp greater than the last imu
% time, I think fimu should be return since it is changing over time
% imuFileType, 0 for text 3dm gx3-35 data, 1 for H764G 1 comprehensive csv data
if(nargin<3)
    imuFileType=0;
end
curimutime=preimutime;
while(preimutime>=curimutime)
    h= fgetl(fimu);
    if (~ischar(h))
        imudata=[];
        return;
    end
    switch(imuFileType)
        case 0
            mass=textscan(h,'%f','delimiter',' ');
            imudata=mass{1};
        case 1
            mass=textscan(h,'%f','delimiter',',');
            imudata=mass{1};
            if(size(imudata,1)<8)
                imudata=[];
                return;
            end
            imudata=imudata(2:8,1);% gps time, xyz delta v, xyz delta theta
            imudata(2:4,1)=imudata(2:4,1)*.3048;% convert to metric unit meter
        case 2
            mass=textscan(h,'%f','delimiter',' ');
            imudata=mass{1};
            hstream= fgetl(fimu);% dummy line
            imu_scalefactor = 9.78;
            imudata(2:4,1) = imudata(2:4,1)/1000.0 * imu_scalefactor;
            imudata(5:7,1) = imudata(5:7,1)*pi/180;
            imudata(6,1) = -imudata(6,1); % wrong gyro y sign in the output
        case 3 % iNEMO output tsv
            mass=textscan(h,'%f','delimiter',',');
            imudata=mass{1};
            imu_scalefactor = 9.8;
            imudata(1)=imudata(1)/1000;
            imudata(2:4,1) = imudata(2:4,1)/1000.0 * imu_scalefactor;
            imudata(5:7,1) = imudata(5:7,1)*pi/180;
            imudata(6,1) = -imudata(6,1); % wrong gyro y sign in the output
        case 4 % microstrain 3dm-gx3-35 csv data file
            mass=textscan(h,'%f','delimiter',',');
            imudata=mass{1};
            imudata=imudata([3,8:13]);
            while(sum(isnan(imudata)))
                h= fgetl(fimu);% GPS line
                if (~ischar(h))
                    imudata=[];
                    return;
                end
                mass=textscan(h,'%f','delimiter',',');
                imudata=mass{1};
                imudata=imudata([3,8:13]);
            end
            imudata(2:4)=imudata(2:4)*9.80665;
        case 5 % epson csv data
            mass=textscan(h,'%f','delimiter',',');
            imudata=mass{1};
            imudata=imudata([2,[8:10, 5:7]]);
            imudata(2:4)=imudata(2:4)*9.80665/1000;
            imudata(5:7)=imudata(5:7)*pi/180;
        otherwise
            fprintf('IMU file type %d is unsupported!\n', imuFileType);
            imudata = [];
            return;
    end
    curimutime=imudata(1,end);
end