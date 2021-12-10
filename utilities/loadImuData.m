function allImuData=loadImuData(imufile, imufiletype, range)
% load imu data in imufile
% output N x M, row format
% time(sec), ax(m/s^2), ay, az, gx (rad/s), gy, gz, and others
if nargin < 3
    range = [0, inf];
end
num_entries=0;
import java.util.LinkedList
preimudata=LinkedList();% record the previous imu data
[fimu, imudata, ~]=readimuheader(imufile, preimudata, range(1), 0, imufiletype); 
allImuData = zeros(10000, 12);
num_entries=num_entries+1;
allImuData(num_entries, 1:length(imudata)) = imudata;
while 1
    [fimu, imudata]=grabnextimudata(fimu, imudata(1), imufiletype);
    if isempty(imudata) || imudata(1) > range(2)
        break;
    else
        num_entries=num_entries+1;
        allImuData(num_entries, 1:length(imudata)) = imudata';
    end
end
allImuData = allImuData(1:num_entries, :);
fclose(fimu);
end
