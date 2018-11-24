function v_rot=rpy(m_rot)
% function v_rot=rpy(m_rot)
%
% Description
%  Computes the orientation vector (in rads) corresponding
%  to a rotation matrix. (dim 3)
%  [psi, theta, phi]
%
% The input datum is:
% - m_rot.- a rotation matrix
%
% The return value is:
%  the rpy vector

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%               Grupo de Robotica                                          %%
%%               Departamento de Informatica e Ingenieria de Sistemas       %%
%%               Universidad de Zaragoza                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Name         : rpy.m                                              %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Programmer   : Jose Maria Martinez                                       %%
%% Languaje     : Matlab 4.2 c                                              %%
%% Date         : February 1996                                             %%
%% Status       : Prueba                                                    %%
%% History	: 16-1-96 creation                                          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



v_rot(3)=angle(m_rot(2,1)*i+m_rot(1,1));
v_rot(2)=angle(-m_rot(3,1)*i+m_rot(1,1)*cos(v_rot(3))+m_rot(2,1)*sin(v_rot(3)));
v_rot(1)=angle(m_rot(3,2)*i+m_rot(3,3));
return;
