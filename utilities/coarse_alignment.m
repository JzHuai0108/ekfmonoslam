function navdata = coarse_alignment(imudata, blh, filename, level, method);
%-------------------------------------------------------
% navdata = coarse_alignment(imudata, gpspos, filename, level);
% Yudan Yi, May 26, 2005
%-------------------------------------------------------
if (exist('level','var')==0 | nargin<4) level = 10; method = 1; end; 
if (isempty(level)) level = 10; end;
if (isempty(method)) method = 2; end;
sample_interval = mean(diff(imudata(:,1)));
sample_rate = floor(1.0/sample_interval+0.5);
if (~isempty(blh))
    if (abs(blh(1))>2*pi && abs(blh(2))>2*pi)   % deg
        blh(1) = blh(1)*pi/180;
        blh(2) = blh(2)*pi/180;
    end
else
    method = 1;
    blh(1) = 39*pi/180;
    blh(2) =-83*pi/180;
    blh(3) = 0.0;
end
[n,m] = size(imudata);
if (level>0)
    for i=2:m
        imudata(:,i) = imuwavedenoise(imudata(:,i), level);
    end
end
index = 0;
numofsec = floor(n/sample_rate);
grav_mag = sqrt(imudata(1,2)*imudata(1,2)+imudata(1,3)*imudata(1,3)+imudata(1,4)*imudata(1,4));
if (grav_mag<0.5*9.8) imudata(:,2:4) = imudata(:,2:4)*sample_rate; end;
switch (method)
    case 1
        att = coarsealign1(imudata,         sample_rate);
    case 2
        att = coarsealign2(imudata, blh(1), sample_rate);
    case 3      
        att = coarsealign3(imudata, blh(1), sample_rate);
    otherwise
        att = coarsealign1(imudata,        sample_rate);
end        
navdata = [att(:,1) ones(size(att,1),1)*blh(1) ones(size(att,1),1)*blh(2) ones(size(att,1),1)*blh(3) zeros(size(att,1),3) att(:,2) att(:,3) att(:,4)];
m = zeros(3,1);
s = zeros(3,1);
for i=1:3
    [m(i),s(i)] = statsangle(navdata(:,i+7));
end
m = m*180/pi;
s = s*180/pi;
% fprintf('%15.10f\t%15.10f\n',[m(:) s(:)]');
% figure
% plot(navdata(:,8:10)*180/pi,'-','Marker','.');
% legend(['\mu=' num2str(m(1),'%5.2f') '\circ \sigma=' num2str(s(1),'%5.2f') '\circ r'],...
%        ['\mu=' num2str(m(2),'%5.2f') '\circ \sigma=' num2str(s(2),'%5.2f') '\circ p'],...
%        ['\mu=' num2str(m(3),'%5.2f') '\circ \sigma=' num2str(s(3),'%5.2f') '\circ y']);
% grid
% xlabel('time');
% ylabel('\circ');
% title('initial coarse alignment for attitude (r,p,h): deg(\circ)');
%     navdata(:, 8) = m(1)*pi/180;
%     navdata(:, 9) = m(2)*pi/180;
%     navdata(:,10) = m(3)*pi/180;
if (~isempty(filename)) 
    fout = fopen(filename,'wb');
    fwrite(fout,navdata','double');
    fclose(fout);
end
