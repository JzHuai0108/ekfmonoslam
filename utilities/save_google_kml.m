function save_google_kml(blh, fileName)
%-------------------------------------------------------------------------------
% save_google_kml(blh, fileName);
% save plots to google map kml format
% input: 
% 1) blh (lat, lon, ht)
% 2) fileName => file name to be saved
% by Dr. Yudan Yi
% modified date on: 04/08/2013
%-------------------------------------------------------------------------------
if (nargin<2) return; end;
[n, m] = size(blh);
if (n<1||m<2) return; end;
fKML = fopen(fileName, 'w');
fprintf(fKML, '<?xml version="1.0" encoding="UTF-8"?>\n');
fprintf(fKML, '<kml xmlns="http://www.opengis.net/kml/2.2">\n');
fprintf(fKML, '<Document>\n');
% fprintf(fKML, '<Placemark>\n');    
% fprintf(fKML, '<name>extruded</name>\n');
% fprintf(fKML, '<LineString>\n');
% fprintf(fKML, '<extrude>1</extrude>\n');
% fprintf(fKML, '<tessellate>1</tessellate>\n');
% fprintf(fKML, '<altitudeMode>relativeToGround</altitudeMode>\n');
% fprintf(fKML, '<coordinates>\n'); 
fprintf(fKML, '<Style id="rov_site">\n');
fprintf(fKML, '<IconStyle>\n');
fprintf(fKML, '<color>ff00f00f</color>\n');
fprintf(fKML, '<scale>0.500</scale>\n');
fprintf(fKML, '<Icon>\n');
fprintf(fKML, '<href>http://maps.google.com/mapfiles/kml/shapes/placemark_circle.png</href>\n');
fprintf(fKML, '</Icon>\n');
fprintf(fKML, '</IconStyle>\n');
fprintf(fKML, '</Style>\n');
for i=1:n
    lat = blh(i,1)*180/pi;
    lon = blh(i,2)*180/pi;
    if m>=3
        ht  = blh(i,3);
    else
        ht = 0.0;
    end
    fprintf(fKML, '<Placemark>\n');
    fprintf(fKML, '<styleUrl>#rov_site</styleUrl>\n');
    fprintf(fKML, '<Point>\n');
    fprintf(fKML, '<coordinates>%14.9f,%14.9f,%14.4f</coordinates>\n', lon, lat, ht);
    fprintf(fKML, '</Point>\n');
    fprintf(fKML, '</Placemark>\n');
%     fprintf(fKML, '%14.9f,%14.9f,%10.3f\n', lon, lat, ht);    
end
% fprintf(fKML, '</coordinates>\n');    
% fprintf(fKML, '</LineString>\n');
% fprintf(fKML, '</Placemark>\n');
fprintf(fKML, '</Document>\n');
fprintf(fKML, '</kml>\n');
fclose(fKML);
