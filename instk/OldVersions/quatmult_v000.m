%quaternian multiplication
%result defines the new rotation
%qab<->Cab, qca<->Cca  => qab*qca=Cab*Cca=Ccb
function qcb=quatmult(qab,qca)
% old implementation
% qcb=zeros(4,1);
% qcb(1)=(qab(1))*qca(1)+(-qab(2))*qca(2)+(-qab(3))*qca(3)+(-qab(4))*qca(4);
% qcb(2)=(qab(2))*qca(1)+(qab(1))*qca(2)+(-qab(4))*qca(3)+(qab(3))*qca(4);
% qcb(3)=(qab(3))*qca(1)+(qab(4))*qca(2)+(qab(1))*qca(3)+(-qab(2))*qca(4);
% qcb(4)=(qab(4))*qca(1)+(-qab(3))*qca(2)+(qab(2))*qca(3)+(qab(1))*qca(4);

qcbx = quaternion(qab) * quaternion(qca);
qcb = compact(qcbx)';
