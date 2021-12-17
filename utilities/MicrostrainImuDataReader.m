
classdef MicrostrainImuDataReader < DataReader
properties(Access=private)
    prevImuBuffer = {};% record the previous imu data 
    imudata = zeros(7, 1);  % time(sec), ax, ay, ax, gx, gy, gz in metric units.
    numPrevImuDataToKeep = 40;
    fileID;
    imuFileType = 4;
    nominalSamplingInterval = 0.008;
end

methods(Access=public)
    function reader = MicrostrainImuDataReader(imufile, startTime)
        [reader.fileID, reader.imudata, reader.prevImuBuffer]=...
            readimuheader(imufile, reader.prevImuBuffer, startTime, ...
            reader.numPrevImuDataToKeep, reader.imuFileType); 
    end
end

methods(Access=protected)
    function data = nextImpl(this)
        if(length(this.prevImuBuffer)==this.numPrevImuDataToKeep) 
            this.prevImuBuffer(1) = [];
        end
        this.prevImuBuffer{end + 1} = this.imudata;
        lasttime = this.imudata(1);
        [this.fileID, this.imudata]=grabnextimudata(...
            this.fileID, this.imudata(1), this.imuFileType);
        currenttime = this.imudata(1);
        if currenttime - lasttime > this.nominalSamplingInterval * 1.5
            fprintf('Warn: Microstrain IMU data missing from %.6f to %.6f!\n', ...
                lasttime, currenttime);
        end
        if currenttime - lasttime < this.nominalSamplingInterval * 0.5
            fprintf('Warn: Microstrain IMU data duplicate from %.6f to %.6f!\n', ...
                lasttime, currenttime);
        end
        data = this.imudata;
    end
    function data = previousImpl(this)
        data = this.prevImuBuffer{end};
    end
    function data = currentImpl(this)
        data = this.imudata;
    end
end
end
