%(Sensor output)=u+Cx+Dv where x_dot=Ax+Bw is the sensor error, u is the
%true kinemativ variable that sensor is expected to measure, Dv is the
%addtive white noise.
%v and w is the WGN with unity spectra
%E{DD'}=R, E{BB'}=Q, sP(sP)'=P0

function [A, B, C, D, sP]=imu_modTI(errdefs)

nsen=length(errdefs);

%form overall system model
A=errdefs(1).A;
B=errdefs(1).B;
C=errdefs(1).C;
D=errdefs(1).D;
sP=errdefs(1).sP;
for in=2:nsen
    A=diagmat_v000(errdefs(in).A,A,[]);
    B=diagmat_v000(errdefs(in).B,B,[]);
    C=diagmat_v000(errdefs(in).C,C,[]);
    D=diagmat_v000(errdefs(in).D,D,[]);
    sP=diagmat_v000(errdefs(in).sP,sP,[]);
end