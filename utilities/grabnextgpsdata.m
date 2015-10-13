function [fgps, gpsdata]=grabnextgpsdata(fgps, gpspostype)
if(~feof(fgps))
    h= fgetl(fgps);
    if (~ischar(h))
        gpsdata(1)=inf;
    else
        if(gpspostype==1)
            mess=sscanf(h,'%d/%d/%d%d:%d:%f%f%f%f%d%d%f%f%f');
            [weeksec, weeknum]=Calender2GPSWeek(mess(1:6));
            gpsdata=[weeksec; mess(7:14)];     
        elseif(gpspostype==2||gpspostype==3)
            stoic=sscanf(h,'%d%f%f%f%f%d%d%f%f%f');
            gpsdata=stoic(2:2+9-1);
        end
    end
    
else
    gpsdata(1)=inf;
end