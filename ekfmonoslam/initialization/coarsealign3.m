function att = coarsealign3(imudata,lat,smoothepoch)
% att = coarsealign3(imudata,lat,smoothepoch);
% initial coarse alignment based on the acceleration and gyro measurement
% given the initial latitude using Dr. Jekeli's book
% Author: Yudan Yi
% Aug. 2004
% Feb. 2005
%--------------------------------------------------------------------------
if (nargin<2) error('error in data'); end;
if (isempty(imudata)|isempty(lat)) error('error in data'); end;
[n,m] = size(imudata);
if (n<1 || m<7) error('error in data'); end;
if (isempty(lat)) error('error in data'); end; 
if (nargin<3) smoothepoch = min(256,n); end;
if (isempty(smoothepoch) || smoothepoch<1) smoothepoch = min(134,n); end;
%--------------------------------------------------------------------------
we = 7.292115147e-05;  % Earth Rotate rate
att = [];
att_index = 0;
smoothepoch=fix(smoothepoch);
for index=1:smoothepoch:n
 	x=index:min(index+smoothepoch-1,n);
    curacc = sum(imudata(x,2:4))/size(imudata,1);
    curgyro = sum(imudata(x,5:7))/size(imudata,1);
    curaccgyro = cross(-curacc,curgyro);
%     gravmag = norm(curaccgyro);
    gravmag = norm(curacc);
    matrix1=[ 0 0 -gravmag; we*cos(lat) 0 -we*sin(lat); 0 gravmag*we*cos(lat) 0];    
    matrix2=[curacc;curgyro;curaccgyro];
    c_nb = matrix1\matrix2;% see Jekeli, 2000, page 243, (8.12)
    curatt = Cbn2att(c_nb);   
    
    att_index = att_index+1;
    att(att_index,:)= [imudata(index,1) curatt(:)'];
end
% [att1, pivot] = extract_pivot_angle(att(:,2:4)');
% meanatt = (mean(att1')'+pivot(:))*180/pi;
% stdatt = (std(att1'))*180/pi;
% fprintf('Mean:\t%15.12f\t%15.12f\t%15.12f\n',meanatt);
% fprintf('Std:\t%15.12f\t%15.12f\t%15.12f\n',stdatt);
%--------------------------------------------------------------------------
% figure
% plot(att(:,2:4)*180/pi,'-','Marker','.');
% legend(['\mu=' num2str(meanatt(1),'%5.2f') '\circ \sigma=' num2str(stdatt(1),'%5.2f') '\circ r'],...
%        ['\mu=' num2str(meanatt(2),'%5.2f') '\circ \sigma=' num2str(stdatt(2),'%5.2f') '\circ p'],...
%        ['\mu=' num2str(meanatt(3),'%5.2f') '\circ \sigma=' num2str(stdatt(3),'%5.2f') '\circ h']);
% grid
% xlabel('time');
% ylabel('\circ');
% title('initial coarse alignment for attitude (h,p,r): deg(\circ)');
% % roll, pitch, heading => heading, pitch and roll
% att(:,2) = meanatt(3)*pi/180; % heading
% att(:,3) = meanatt(2)*pi/180; % pitch
% att(:,4) = meanatt(1)*pi/180; % roll