%dcm (Cba) to quaternion (qba)
%Cba:B to A transformation matrix
%qba:B to A quaternion

%Implementation of savage(3-47)
function quat=dcm2quat(dcm)

tr=dcm(1,1)+dcm(2,2)+dcm(3,3);
Pa=1+tr;
Pb=1+2*dcm(1,1)-tr;
Pc=1+2*dcm(2,2)-tr;
Pd=1+2*dcm(3,3)-tr;

quat=zeros(4,1);
if (Pa>=Pb && Pa>=Pc && Pa>=Pd)
    quat(1)=0.5*sqrt(Pa);
    quat(2)=(dcm(3,2)-dcm(2,3))/4/quat(1);
    quat(3)=(dcm(1,3)-dcm(3,1))/4/quat(1);
    quat(4)=(dcm(2,1)-dcm(1,2))/4/quat(1);
elseif (Pb>=Pc && Pb>=Pd)
    quat(2)=0.5*sqrt(Pb);
    quat(3)=(dcm(2,1)+dcm(1,2))/4/quat(2);
    quat(4)=(dcm(1,3)+dcm(3,1))/4/quat(2);
    quat(1)=(dcm(3,2)-dcm(2,3))/4/quat(2);
elseif (Pc>=Pd)
    quat(3)=0.5*sqrt(Pc);
    quat(4)=(dcm(3,2)+dcm(2,3))/4/quat(3);
    quat(1)=(dcm(1,3)-dcm(3,1))/4/quat(3);
    quat(2)=(dcm(2,1)+dcm(1,2))/4/quat(3);
else
    quat(4)=0.5*sqrt(Pd);
    quat(1)=(dcm(2,1)-dcm(1,2))/4/quat(4);
    quat(2)=(dcm(1,3)+dcm(3,1))/4/quat(4);
    quat(3)=(dcm(3,2)+dcm(2,3))/4/quat(4);
end

if (quat(1)<=0)
    quat=-quat;
end

