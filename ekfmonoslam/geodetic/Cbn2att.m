function att = Cbn2att(Cbn);
%-------------------------------------------------------
% Yudan Yi, May 26, 2005
%-------------------------------------------------------
% attitude: roll, pitch and heading
att(1) = atan2(Cbn(3,2),Cbn(3,3));
att(2) = asin(-Cbn(3,1));
att(3) = atan2(Cbn(2,1),Cbn(1,1));
att = att(:);
return;