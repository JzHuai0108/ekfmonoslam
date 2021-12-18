function Cba=quat2dcm(qba)

a=qba(1);
b=qba(2);
c=qba(3);
d=qba(4);

Cba=zeros(3);
Cba(1,1)=a^2+b^2-c^2-d^2;
Cba(2,1)=2*(b*c+a*d);
Cba(3,1)=2*(b*d-a*c);

Cba(1,2)=2*(b*c-a*d);
Cba(2,2)=a^2-b^2+c^2-d^2;
Cba(3,2)=2*(c*d+a*b);

Cba(1,3)=2*(b*d+a*c);
Cba(2,3)=2*(c*d-a*b);
Cba(3,3)=a^2-b^2-c^2+d^2;

