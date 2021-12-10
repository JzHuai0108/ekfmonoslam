function gpsdata_all=loadAllGPSData(gpsfile, gpsSE, coordinateStyle)
% load rtklib output of gps solutions, input gpsfile name, gps start and
% end time of gps TOW, gpspostype are as described in readgpsheader
% output: Nx9 gps entries in the session gpsSE, each entry is 
% GPS TOW, ECEF XYZ, Q and no of satels and sdx sdy sdz

% if coordinateStyle == lla, the first four columns will be 
% GPS TOW, lat lon ellipsoid height. 
if nargin < 3 
    coordinateStyle = 'ecef'; 
end

gpsdata_all=zeros( ceil((gpsSE(2)- gpsSE(1)+100)*5), 12);
num_entries=0;
[fgps, gpsdata, gpspostype]=readgpsheader(gpsfile, gpsSE(1));
num_entries=num_entries+1;
gpsdata_all(num_entries, :) = gpsdata';
while(1)
    [fgps, gpsdata]=grabnextgpsdata(fgps, gpspostype);
    if(gpsdata(1)==inf || gpsdata(1)> gpsSE(2))
        break;
    else
        num_entries=num_entries+1;
        gpsdata_all(num_entries, :) = gpsdata';
    end
end
fclose(fgps);
gpsdata_all= gpsdata_all(1:num_entries, :);
if strcmp(coordinateStyle, 'lla')
    gpsdata_all(:, 2:4) = ecef2lla(gpsdata_all(:, 2:4));
    gpsdata_all(:, 2:3) = gpsdata_all(:, 2:3) * pi / 180;
end
end