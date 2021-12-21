function att = qua2att(qua);
%-------------------------------------------------------
% att = qua2att(qua);
% Yudan Yi, May 26, 2005
%-------------------------------------------------------
Cbn = qua2Cbn(qua);
att = Cbn2att(Cbn);
att = att(:);
return;
