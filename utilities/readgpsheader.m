function [fgps, gpsdata]=readgpsheader(gpsfile, preimutime, gpspostype)
fgps=fopen(gpsfile,'rt');
if(fgps~=-1)
    %%%%%Discard all observation data before the current time
    %remove the header of gps posdata
    hstream= fgetl(fgps);
    while(true)
        if(isempty(strfind(hstream,'%')))
            break;
        else hstream= fgetl(fgps);
        end
    end
    gpsdata=zeros(1,9); %GPS TOW, lattitude longitude in degree and height in meter, Q and no of satels and sdn sde sdu
    if(gpspostype==1)
        stoic=sscanf(hstream,'%d/%d/%d%d:%d:%f%f%f%f%d%d%f%f%f');
        [weeksec, weeknum]=Calender2GPSWeek(stoic(1:6));
        mess=stoic;
        
        gpsdata(1)=weeksec+(mess(4:6)-stoic(4:6))'*[3600;60;1];
        gpsdata(2:9)=mess(7:14);        
        while (gpsdata(1)<=preimutime)
            hstream= fgetl(fgps);
            mess=sscanf(hstream,'%d/%d/%d%d:%d:%f%f%f%f%d%d%f%f%f');
            gpsdata(1)=weeksec+(mess(4:6)-stoic(4:6))'*[3600;60;1];
            gpsdata(2:9)=mess(7:14);
        end
    elseif(gpspostype==2)
        stoic=sscanf(hstream,'%d%f%f%f%f%d%d%f%f%f');
        gpsdata=stoic(2:2+9-1); % GPS TOW lattitude longitude in degree and height in meter, Q and no of satels and std of N,E,U
        
        while (gpsdata(1)<=preimutime)
            hstream= fgetl(fgps);
            stoic=sscanf(hstream,'%d%f%f%f%f%d%d%f%f%f');
            gpsdata=stoic(2:2+9-1);
        end
    elseif(gpspostype==3)
        stoic=sscanf(hstream,'%d%f%f%f%f%d%d%f%f%f');
        gpsdata=stoic(2:2+9-1); % GPS TOW ecef X,Y,Z, Q and no of satels and std ecef X,Y,Z
        
        while (gpsdata(1)<=preimutime)
            hstream= fgetl(fgps);
            stoic=sscanf(hstream,'%d%f%f%f%f%d%d%f%f%f');
            gpsdata=stoic(2:2+9-1);
        end
    end
else gpsdata=inf;
end