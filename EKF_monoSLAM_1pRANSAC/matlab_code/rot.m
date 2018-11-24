function r=rot(psi,theta,phi);
% function r=rot(psi,theta,phi)
%
% Description
%  Computes the rotarion matrix corresponding to a oritentation vector
% The input data are:
% - psi   .- the psi angle (radians)
% - theta .- the theta angle (radians)
% - phi   .- the phi angle (radians)
%
% The return value is
% the rotation matrix (dim 3x3)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%               Grupo de Robotica                                          %%
%%               Departamento de Informatica e Ingenieria de Sistemas       %%
%%               Universidad de Zaragoza                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Name         : .m                                              %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Programmer   : Jose Maria Martinez                                       %%
%% Languaje     : Matlab 4.2 c                                              %%
%% Date         : February 1996                                             %%
%% Status       : Prueba                                                    %%
%% History	: 16-8-95 creation                                          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





cosphi=cos(phi);
sinphi=sin(phi);
costheta=cos(theta);
sintheta=sin(theta);
cospsi=cos(psi);
sinpsi=sin(psi);

r=zeros(3,3);
r(1,1)=cosphi*costheta;
r(1,2)=cosphi*sintheta*sinpsi-sinphi*cospsi;
r(1,3)=cosphi*sintheta*cospsi+sinphi*sinpsi;
r(2,1)=sinphi*costheta;
r(2,2)=sinphi*sintheta*sinpsi+cosphi*cospsi;
r(2,3)=sinphi*sintheta*cospsi-cosphi*sinpsi;
r(3,1)=-sintheta;
r(3,2)=costheta*sinpsi;
r(3,3)=costheta*cospsi;
return;
