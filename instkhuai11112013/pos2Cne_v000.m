%given lat and lon computes ned to e dcm
function Cne=pos2Cne(lat,lon)

sL=sin(lat);
cL=cos(lat);
sl=sin(lon);
cl=cos(lon);

Cne=[-sL*cl -sl -cL*cl;
     -sL*sl cl  -cL*sl;
     cL     0   -sL];