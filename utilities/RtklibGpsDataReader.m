
classdef RtklibGpsDataReader < DataReader
properties(Access=private)
    gpsdata = zeros(12, 1);  % time(sec), x, y, z, Q, n, cov x, cov y, cov z, and others in metric units.
    fileID;
    gnssFileType = 0;
    nominalSamplingInterval = 0.2;
    outputFormat = 'lla';
end

methods(Access=public)
    function reader = RtklibGpsDataReader(gpsfile, startTime, outputFormat)
        reader.outputFormat = outputFormat;
        [reader.fileID, reader.gpsdata, reader.gnssFileType]=...
            readgpsheader(gpsfile, startTime, outputFormat); 
    end
end

methods(Access=protected)
    function data = nextImpl(this)
        lasttime = this.gpsdata(1);
        [this.fileID, this.gpsdata]=grabnextgpsdata(this.fileID, this.gnssFileType, this.outputFormat);
        currenttime = this.gpsdata(1);
        if currenttime - lasttime > this.nominalSamplingInterval * 1.5
            fprintf('Warn: RTKLib GNSS data missing from %.6f to %.6f!\n', ...
                lasttime, currenttime);
        end
        if currenttime - lasttime < this.nominalSamplingInterval * 0.5
            fprintf('Warn: RTKLib GNSS data duplicate from %.6f to %.6f!\n', ...
                lasttime, currenttime);
        end
        data = this.gpsdata;
    end
    function data = previousImpl(this)
        data = [];
    end
    function data = currentImpl(this)
        data = this.gpsdata;
    end
end
end