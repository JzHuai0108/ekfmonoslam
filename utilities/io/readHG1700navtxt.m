function navdata=readHG1700navtxt(navFileName, deltaT)
% readHG1700imunavtxt navfilename, for HG1700, deltaT data duration 
%in seconds %read in the imunav.txt of HG1700
if (isempty(navFileName)) return; end;
fpos = fopen(navFileName,'r'); if (fpos==-1) return; end;
stk=(deltaT+60)*100;% 60s is enough
navdata=zeros(stk,7);
i=0;
while (i<stk)  
    g = fscanf(fpos,'%e', 7);    
    if (isempty(g))
        break
    end
    h  = fscanf(fpos,'%s', 3); 
    k  = fscanf(fpos,'%e', 6);
    i=i+1;
    navdata(i,:) = [g(1) k'];    
end
fclose(fpos);
navdata=navdata(1:i,:);
navdata(:,5:7)=navdata(:,5:7)*0.3048;
temp=navdata(:,5:7);
navdata(:,5:7)=navdata(:,2:4);
navdata(:,2:4)=temp;