function [fgps, gpsdata]=grabnextgpsdata(fgps, gpspostype, coordinateStyle)
% output: gpsdata 9x1 GPS TOW, ECEF XYZ, Q and no of satels and sdx sdy sdz
% sdxy sdyz sdzx
% the gpspostype variable is determined based on the first line of data 
% (whether it is GPST or GPS TOW) and last line of the header to determine 
% the coordinated type (this line contains unique key words ecef or deg)
% more information about gpspostype in readgpsheader
if nargin < 3
    coordinateStyle = 'ecef';
end
gpsdata = zeros(12, 1);
if ~feof(fgps)
    hstream= fgetl(fgps);
    if (~ischar(hstream))
        gpsdata(1)=inf;
    else
        if(gpspostype==1 || gpspostype==4 || gpspostype==5)
            if (gpspostype==5)% d'"
                mess=sscanf(hstream,'%d/%d/%d%d:%d:%f%d%d%f%d%d%f%f%d%d%f%f%f%f%f%f');
                latdeg=sign(mess(7))*sum(abs(mess(7:9)).*[1;1/60;1/3600]);
                londeg=sign(mess(10))*sum(abs(mess(10:12)).*[1;1/60;1/3600]);
                mess=[mess(1:6);latdeg;londeg;mess(13:end)];
            else
                mess=sscanf(hstream,'%d/%d/%d%d:%d:%f%f%f%f%d%d%f%f%f%f%f%f');
            end
            [weeksec, ~]=Calender2GPSWeek(mess(1:6));
            gpsdata=[weeksec; mess(7:17)];
            if(gpspostype==1 || gpspostype==5)
                %convert degree llh to ecef xyz
                gpsdata(2:4)= ecef2geo_v000([gpsdata(2)/180*pi;gpsdata(3)/180*pi;gpsdata(4)],1)';
                % convert sdn,e,u to sdx,y,z
                Ce2n=llh2dcm_v000([gpsdata(2)/180*pi;gpsdata(3)/180*pi;gpsdata(4)]);
                covm=cov2RTKlib(gpsdata(7:12),1); % covariance matrix from RTKlib notation
                gpsdata(7:12)=cov2RTKlib(Ce2n'*covm*Ce2n,0); % back to RTKlib notation but transformed
            end
        elseif(gpspostype==2||gpspostype==3 || gpspostype==6)
            if (gpspostype==6)% d'"
                stoic=sscanf(hstream,'%d%f%d%d%f%d%d%f%f%d%d%f%f%f%f%f%f');
                latdeg=sign(stoic(3))*sum(abs(stoic(3:5)).*[1;1/60;1/3600]);
                londeg=sign(stoic(6))*sum(abs(stoic(6:8)).*[1;1/60;1/3600]);
                stoic=[stoic(1:2);latdeg;londeg;stoic(9:end)];
            else
                stoic=sscanf(hstream,'%d%f%f%f%f%d%d%f%f%f%f%f%f');
            end
            gpsdata=stoic(2:2+12-1);
            if(gpspostype==2 || gpspostype==6)
                %convert degree llh to ecef xyz
                gpsdata(2:4)= ecef2geo_v000([gpsdata(2)/180*pi;gpsdata(3)/180*pi;gpsdata(4)],1)';
                % convert sdn,e,u to sdx,y,z
                Ce2n=llh2dcm_v000([gpsdata(2)/180*pi;gpsdata(3)/180*pi;gpsdata(4)]);
                covm=cov2RTKlib(gpsdata(7:12),1); % covariance matrix from RTKlib notation
                gpsdata(7:12)=cov2RTKlib(Ce2n'*covm*Ce2n,0); % back to RTKlib notation but transformed
            end
        else
            fprintf('GNSS file type %d is unsupported!\n', gpspostype);
            gpsdata(1) = inf;
        end
    end
else
    gpsdata(1)=inf;
end
if strcmp(coordinateStyle, 'lla') && gpsdata(1) ~= inf
    gpsdata(2:4) = ecef2lla(gpsdata(2:4)')';
    gpsdata(2:3) = gpsdata(2:3) * pi / 180;
end
end