%dir=0 -> q=qab (vb=Cab*va)
%dir=1 -> q=qba (vb=Cba'*va)
function vb=quatrot(qab,va,dir)

va=[0;va(1);va(2);va(3)];
if dir==0
    q=qab;
else
    q=qab;
    q(2:4)=-q(2:4);
end
vr_a=quatmult_v001(q,va,0);
q(2:4)=-q(2:4);
vb=quatmult_v001(vr_a,q,0);
vb=vb(2:4);
