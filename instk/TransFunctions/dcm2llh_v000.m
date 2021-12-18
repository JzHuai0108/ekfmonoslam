%wander=[sin(alpha);cos(alpha)]
function [llh, wander]=dcm2llh(Cen)

llh=zeros(2,1); %despite the name, no height
llh(1)=-asin(Cen(3,3));
cL=cos(llh(1));
llh(2)=atan2(-Cen(3,2)/cL, -Cen(3,1)/cL);

if (nargout==2)
    wander=zeros(2,1);
    wander(1)=-Cen(2,3)/cL;
    wander(2)=Cen(1,3)/cL;
end
