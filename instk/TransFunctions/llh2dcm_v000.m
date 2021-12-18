%wander=[sin(alpha);cos(alpha)] (optional)
%(rotate geo frame around z with "wander angle" to obtain nav frame)
function Cen=llh2dcm(llh, wander)

if nargin==1
    wander=[0;1];
end

if (isempty(wander))
    wander=zeros(2,1);
    wander(2)=1;
end

sL=sin(llh(1));
cL=cos(llh(1));
sl=sin(llh(2));
cl=cos(llh(2));

Cgn=[wander(2) wander(1) 0;-wander(1) wander(2) 0;0 0 1];
Ceg=[-sL*cl -sL*sl cL;-sl cl 0;-cL*cl -cL*sl -sL];

Cen=Cgn*Ceg;