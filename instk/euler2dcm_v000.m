%eul defined in "n": rotate "n" to obtain "b"
%result: Cbn (from b to n)
function dcm = euler2dcm(eul)
cr=cos(eul(1));sr=sin(eul(1));	%roll
cp=cos(eul(2));sp=sin(eul(2));	%pitch
ch=cos(eul(3));sh=sin(eul(3));	%heading
dcm = zeros(3);

dcm(1,1)=cp*ch;
dcm(1,2)=(sp*sr*ch)-(cr*sh);
dcm(1,3)=(cr*sp*ch)+(sh*sr);

dcm(2,1)=cp*sh;
dcm(2,2)=(sr*sp*sh)+(cr*ch);
dcm(2,3)=(cr*sp*sh)-(sr*ch);

dcm(3,1)=-sp;
dcm(3,2)=sr*cp;
dcm(3,3)=cr*cp;