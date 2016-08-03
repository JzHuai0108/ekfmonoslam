function gpsdata_all=loadAllGPSData(gpsfile, gpsSE, gpspostype)
% load rtklib output of gps solutions, input gpsfile name, gps start and
% end time of gps TOW, gpspostype are as described in readgpsheader
% output: Nx9 gps entries in the session gpsSE, each entry is 
% GPS TOW, ECEF XYZ, Q and no of satels and sdx sdy sdz

gpsdata_all=zeros( ceil((gpsSE(2)- gpsSE(1)+100)*5), 12);
num_entries=0;
[fgps, gpsdata]=readgpsheader(gpsfile, gpsSE(1)-0.2, gpspostype);
num_entries=num_entries+1;
gpsdata_all(num_entries, :) = gpsdata';
while(1)
    [fgps, gpsdata]=grabnextgpsdata(fgps, gpspostype);
    if(gpsdata(1)==inf || gpsdata(1)> gpsSE(2)+0.2)
        break;
    else
        num_entries=num_entries+1;
        gpsdata_all(num_entries, :) = gpsdata';
    end
end
gpsdata_all= gpsdata_all(1:num_entries, :);
fclose(fgps);
end