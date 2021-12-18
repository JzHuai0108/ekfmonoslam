function qbn = att2qua(att);
%-------------------------------------------------------
% qua = att2qua(att);
% Yudan Yi, May 26, 2005
% Huai. 10/2/2013
%-------------------------------------------------------
Cbn=euler2dcm_v000(att);
qbn = dcm2quat(Cbn);
