%if v==1 --> res=qab'*qca
%if v==2 --> res=qab*qca'
%else res=qab*qca=qcb
function res=quatmult(qab,qca, v)

% res=zeros(4,1);
% if (v==1)
%     qab=[qab(1);-qab(2);-qab(3);-qab(4)];
% elseif (v==2)
%     qca=[qca(1);-qca(2);-qca(3);-qca(4)];
% end
% 
% res(1)=(qab(1))*qca(1)+(-qab(2))*qca(2)+(-qab(3))*qca(3)+(-qab(4))*qca(4);
% res(2)=(qab(2))*qca(1)+(qab(1))*qca(2)+(-qab(4))*qca(3)+(qab(3))*qca(4);
% res(3)=(qab(3))*qca(1)+(qab(4))*qca(2)+(qab(1))*qca(3)+(-qab(2))*qca(4);
% res(4)=(qab(4))*qca(1)+(-qab(3))*qca(2)+(qab(2))*qca(3)+(qab(1))*qca(4);
% 
if v == 1
    return compact(conj(quaternion(qab)) * quaternion(qca));
elseif v == 2
    return compact(quaternion(qab) * conj(quaternion(qca)));
else
    return compact(quaternion(qab) * quaternion(qca));
end
