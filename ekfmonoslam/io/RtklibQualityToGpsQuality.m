function gpsQuality = RtklibQualityToGpsQuality(rtklibquality)
gpsQualityList = [4,5,9,2,1,5];
if rtklibquality < 7
    gpsQuality = gpsQualityList(rtklibquality);
else
    gpsQuality = 0;
end
end
