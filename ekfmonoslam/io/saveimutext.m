function saveimutext(filename, table, header)
%to save the data in the table with a header string to the file of iflename
%flush all previous data in the file if it exists
fileID = fopen(filename,'w');
header='%GPS TOW, AX, AY, AZ, WX, WY, WZ';
fprintf(fileID,'%s\n',header);
[m,n]=size(table);
for i=1:m
fprintf(fileID, '%15.9f,%18.15g,%18.15g,%18.15g,%18.15g,%18.15g,%18.15g\n', table(i,:)); 
% save with 15 significant digits with %g,  '%6.4g' prints pi as ' 3.142'
end
fclose(fileID);