% rovXYZ 3 x N vectors, inillh 3 x N vectors
% output 3 x N vectors
function dNED=posdiff_v001(rovXYZ, inillh)
% for blh2xyz, input N x 3, output N x 3;
startXYZ= blh2xyz(inillh');
startCen=llh2dcm_v000(inillh(1:2),[0,1]);
dXYZ = rovXYZ-startXYZ';
dNED = startCen*dXYZ;