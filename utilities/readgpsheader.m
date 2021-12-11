function [fgps, gpsdata, gpspostype]=readgpsheader(gpsfile, preimutime, coordinateStyle)
% return the first gpsdata as a column vector that has a timestamp
% greater than preimutime, and fgps the file Pointer to the gpsfile
% gpspostype variable determines the type of the RTKlib file:
% 1 GPS h/m/s + lat lon h (deg deg m)
% 2 GPS TOW + lat lon h (deg deg m)
% 3 GPS TOW + ECEF XYZ
% 4 GPS h/m/s + ECEF XYZ
% 5 GPST h/m/s, lat lon h (d ' " d ' " m)
% 6 GPS TOW, lat lon h (d ' " d ' " m)
% coordinateStyle == 'ecef' or 'lla'

% output: gpsdata 12x1 GPS TOW, ECEF XYZ (or lla depending on coordinateStyle), 
% Q and no of satels, sdx sdy sdz,
% sdxy, sdyz, sdzx (according the RTKlib scheme; conversion to 3x3 matrix uses covm2RTKlib)
if nargin < 3
    coordinateStyle = 'ecef';
end
gpsdata=zeros(12, 1);
fgps=fopen(gpsfile,'rt');
if fgps~=-1
    %% Discard all observation data before the current time
    %remove the header of gps posdata
    hstream = fgetl(fgps);
    hstream0 = hstream;
    while(true)
        if ~contains(hstream,'%')
            break;
        else
            hstream0=hstream; % keeping old line to recognize type
            hstream= fgetl(fgps);
        end
    end
    
    %% Recognizing the type of the RTKlib pos file and owerwite the gpspostype variable
    % checking the time type
    if ~contains(hstream,'/')
        time_type=1; % GPS TOW
    else
        time_type=0; % GPST
    end
    % checking the data type
    if ~contains(hstream0,'ecef')
        if ~contains(hstream0,'deg')
            coord_type=5; % d'" for lat and lon
        else
            coord_type=1; % deg for lat and lon
        end
    else
         coord_type=2; % ecef coords
    end
    gpspostype=time_type+coord_type; % should match previously written values
    if time_type==0 && coord_type==2 % one needs to be override
        gpspostype=4;
    end
    
    %% reading data
    if(gpspostype==1 || gpspostype==5)
        if (gpspostype==5)% d'"
            stoic=sscanf(hstream,'%d/%d/%d%d:%d:%f%d%d%f%d%d%f%f%d%d%f%f%f%f%f%f');
            latdeg=sign(stoic(7))*sum(abs(stoic(7:9)).*[1;1/60;1/3600]);
            londeg=sign(stoic(10))*sum(abs(stoic(10:12)).*[1;1/60;1/3600]);
            stoic=[stoic(1:6);latdeg;londeg;stoic(13:end)];
        else
            stoic=sscanf(hstream,'%d/%d/%d%d:%d:%f%f%f%f%d%d%f%f%f%f%f%f');
        end
        %GPST, lattitude longitude in degree and height in meter, Q and no of satels and sdn sde sdu
        [weeksec, ~]=Calender2GPSWeek(stoic(1:6));
        mess=stoic;
        
        gpsdata(1)=weeksec+(mess(4:6)-stoic(4:6))'*[3600;60;1];
        gpsdata(2:12)=mess(7:17);
        %convert degree llh to ecef xyz
        gpsdata(2:4)= ecef2geo_v000([gpsdata(2)/180*pi;gpsdata(3)/180*pi;gpsdata(4)],1)';
        % convert sdn,e,u to sdx,y,z
        Ce2n=llh2dcm_v000([gpsdata(2)/180*pi;gpsdata(3)/180*pi;gpsdata(4)]);
        covm=cov2RTKlib(gpsdata(7:12),1); % covariance matrix from RTKlib notation
        gpsdata(7:12)=cov2RTKlib(Ce2n'*covm*Ce2n,0); % back to RTKlib notation but transformed
        while (gpsdata(1)<=preimutime)
            hstream= fgetl(fgps);
            if (gpspostype==5)% d'"
                mess=sscanf(hstream,'%d/%d/%d%d:%d:%f%d%d%f%d%d%f%f%d%d%f%f%f%f%f%f');
                latdeg=sign(mess(7))*sum(abs(mess(7:9)).*[1;1/60;1/3600]);
                londeg=sign(mess(10))*sum(abs(mess(10:12)).*[1;1/60;1/3600]);
                mess=[mess(1:6);latdeg;londeg;mess(13:end)];
            else
                mess=sscanf(hstream,'%d/%d/%d%d:%d:%f%f%f%f%d%d%f%f%f%f%f%f');
            end
            gpsdata(1)=weeksec+(mess(4:6)-stoic(4:6))'*[3600;60;1];
            gpsdata(2:12)=mess(7:17);
            %convert degree llh to ecef xyz
            gpsdata(2:4)= ecef2geo_v000([gpsdata(2)/180*pi;gpsdata(3)/180*pi;gpsdata(4)],1)';
            % convert sdn,e,u to sdx,y,z
            Ce2n=llh2dcm_v000([gpsdata(2)/180*pi;gpsdata(3)/180*pi;gpsdata(4)]);
            covm=cov2RTKlib(gpsdata(7:12),1); % covariance matrix from RTKlib notation
            gpsdata(7:12)=cov2RTKlib(Ce2n'*covm*Ce2n,0); % back to RTKlib notation but transformed
        end
    elseif(gpspostype==2 || gpspostype==6)
        if (gpspostype==6)% d'"
            stoic=sscanf(hstream,'%d%f%d%d%f%d%d%f%f%d%d%f%f%f%f%f%f');
            latdeg=sign(stoic(3))*sum(abs(stoic(3:5)).*[1;1/60;1/3600]);
            londeg=sign(stoic(6))*sum(abs(stoic(6:8)).*[1;1/60;1/3600]);
            stoic=[stoic(1:2);latdeg;londeg;stoic(9:end)];
        else
            stoic=sscanf(hstream,'%d%f%f%f%f%d%d%f%f%f%f%f%f');
        end
        % GPS TOW lattitude longitude in degree and height in meter, Q and no of satels and std of N,E,U
        gpsdata=stoic(2:2+12-1);
        %convert degree llh to ecef xyz
        gpsdata(2:4)= ecef2geo_v000([gpsdata(2)/180*pi;gpsdata(3)/180*pi;gpsdata(4)],1)';
        % convert sdn,e,u to sdx,y,z
        Ce2n=llh2dcm_v000([gpsdata(2)/180*pi;gpsdata(3)/180*pi;gpsdata(4)]);
        covm=cov2RTKlib(gpsdata(7:12),1); % covariance matrix from RTKlib notation
        gpsdata(7:12)=cov2RTKlib(Ce2n'*covm*Ce2n,0); % back to RTKlib notation but transformed
        while (gpsdata(1)<=preimutime)
            hstream= fgetl(fgps);
            if (gpspostype==6)% d'"
                stoic=sscanf(hstream,'%d%f%d%d%f%d%d%f%f%d%d%f%f%f%f%f%f');
                latdeg=sign(stoic(3))*sum(abs(stoic(3:5)).*[1;1/60;1/3600]);
                londeg=sign(stoic(6))*sum(abs(stoic(6:8)).*[1;1/60;1/3600]);
                stoic=[stoic(1:2);latdeg;londeg;stoic(9:end)];
            else
                stoic=sscanf(hstream,'%d%f%f%f%f%d%d%f%f%f%f%f%f');
            end
            gpsdata=stoic(2:2+12-1);
            %convert degree llh to ecef xyz
            gpsdata(2:4)= ecef2geo_v000([gpsdata(2)/180*pi;gpsdata(3)/180*pi;gpsdata(4)],1)';
            % convert sdn,e,u to sdx,y,z
            Ce2n=llh2dcm_v000([gpsdata(2)/180*pi;gpsdata(3)/180*pi;gpsdata(4)]);
            covm=cov2RTKlib(gpsdata(7:12),1); % covariance matrix from RTKlib notation
            gpsdata(7:12)=cov2RTKlib(Ce2n'*covm*Ce2n,0); % back to RTKlib notation but transformed
        end
    elseif(gpspostype==3)
        stoic=sscanf(hstream,'%d%f%f%f%f%d%d%f%f%f%f%f%f');% GPS TOW ecef X,Y,Z, Q and no of satels and std ecef X,Y,Z
        gpsdata=stoic(2:2+12-1);
        
        while (gpsdata(1)<=preimutime)
            hstream= fgetl(fgps);
            stoic=sscanf(hstream,'%d%f%f%f%f%d%d%f%f%f%f%f%f');
            gpsdata=stoic(2:2+12-1);
        end
       
    elseif(gpspostype==4)
        %  GPST    x-ecef(m)      y-ecef(m)      z-ecef(m)   Q  ns
        % sdx(m)   sdy(m)   sdz(m)  sdxy(m)  sdyz(m)  sdzx(m) age(s)  ratio
        stoic=sscanf(hstream,'%d/%d/%d%d:%d:%f%f%f%f%d%d%f%f%f%f%f%f');
        [weeksec, weeknum]=Calender2GPSWeek(stoic(1:6));
        mess=stoic;
        %gpsdata is GPS TOW ecef X,Y,Z, Q and no of satels and std ecef X,Y,Z
        gpsdata(1)=weeksec+(mess(4:6)-stoic(4:6))'*[3600;60;1];
        gpsdata(2:12)=mess(7:17);
        while (gpsdata(1)<=preimutime)
            hstream= fgetl(fgps);
            mess=sscanf(hstream,'%d/%d/%d%d:%d:%f%f%f%f%d%d%f%f%f%f%f%f');
            gpsdata(1)=weeksec+(mess(4:6)-stoic(4:6))'*[3600;60;1];
            gpsdata(2:12)=mess(7:17);
        end
    else
        fprintf('GNSS file type %d is unsupported!\n', gpspostype);
    end
else
    fprintf('Failed to open GNSS file at %s.\n', gpsfile);
    gpsdata(1)=inf;
end
if strcmp(coordinateStyle, 'lla') && gpsdata(1) ~= inf
    gpsdata(2:4) = ecef2lla(gpsdata(2:4)')';
    gpsdata(2:3) = gpsdata(2:3) * pi / 180;
end
end

function testReadGPSHeader()
    filename= 'E:\Novatel_rover.pos';
    [filePointer, gpsData]=readgpsheader(filename,  316799.2, 4)
    [filePointer, gpsData]=grabnextgpsdata(filePointer, 4)
    fclose(filePointer);
    filename2 ='F:\20130808\oem615_20130809.pos';
    [filePointer, gpsData]=readgpsheader(filename2,   414790.4, 2)
    [filePointer, gpsData]=grabnextgpsdata(filePointer, 2)
    fclose(filePointer);
    filename3 ='F:\3_Nikon\GPS_OEM615-Nikon\Processed_respect_Van_base_station\Nikon.pos';
    [filePointer, gpsData]=readgpsheader(filename3, 324595.4, 3)
    [filePointer, gpsData]=grabnextgpsdata(filePointer, 3)
    fclose(filePointer);
end