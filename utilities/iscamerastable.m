% statehistory is of CQueue()
% radtheshold the threshold for angles in unit of radians
% distthreshold the threshold for distance changes in unit of meters
function [judge, estimate]=iscamerastable(databuf, radtheshold, distthreshold)
rows=size(databuf,1);
cols=size(databuf,2);
estimate=zeros(rows,1);
% obtain the average for quaternions
estimate(1:4)=weighted_avg_quat(databuf(1:4,:)',ones(cols,1));
estimate(5:7)=mean(databuf(5:7,:),2);
delta=zeros(6,cols);
for inst=1:cols
delta(4:6, inst)=databuf(5:7, inst)-estimate(5:7);
% kind of like subangle(qua2att(X1(7:10,i)),qua2att(X2(7:10,i)))
delta(1:3, inst)=quat2rot_v000(quatmult_v001(databuf(1:4,inst),estimate(1:4),2));
end
stdevs=sqrt(sum(delta.*delta,2)/cols);
% later use Mahalanobis distance would be better
judge=false;
if(max(stdevs(1:3))<radtheshold&&max(stdevs(4:6))<distthreshold)
    judge=true;
end