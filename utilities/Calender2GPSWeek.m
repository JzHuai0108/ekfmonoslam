function [weekSec, weekNum]=Calender2GPSWeek(bundle)

year=bundle(1);
month=bundle(2);
day=bundle(3);
hour=bundle(4);
minute=bundle(5);
second=bundle(6);
if (year < 80)
    year =year+ 2000;
elseif (year > 80 && year < 1900)
    year =year+ 1900;
end
totalDay=0;
if (year>=1981)
    totalDay=360;
end
for yearIndex=1981:1:(year-1)
    totalDay =totalDay+ 365;
    if ( (mod(yearIndex,4)==0&&mod(yearIndex,100)~=0)||mod(yearIndex,400)==0 )
        totalDay=totalDay+1;
    end
end
dayPerMon=[31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
for ( monIndex=1:1:(month-1))
    totalDay =totalDay+ dayPerMon(monIndex);
end
if (month>2 && ((mod(year,4)==0&&mod(year,100)~=0)||mod(year,400)==0))
    totalDay=totalDay+1;
end
totalDay=totalDay+day;
weekNum	=floor(totalDay/7);
weekSec	= (totalDay-weekNum*7)*24.0*3600.0+hour*3600.0+minute*60.0+second;



