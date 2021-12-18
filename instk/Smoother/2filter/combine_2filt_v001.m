%Combine forward and backward solutions

%Note that the difference between this and the combine_2filt_v000 is very
%minimal. Therefore, you can use any of these. 

%dzBck=PbckInv*dx_bckwrd
function [P_smt, Cen_smt, h_smt, Vn_smt, qbn_smt, imu_smt, dxSmt]=combine_2filt(fwdnav, fwdcov, Cen_bck, h_bck, Vn_bck, qbn_bck, dzBck, PbckInv)

%Forward filter results
Cen_fwd=reshape(fwdnav(1:9),[3,3]);
h_fwd=fwdnav(10); %elipsoid height
Vn_fwd=fwdnav(11:13);
qbn_fwd=fwdnav(14:17);
imu_fwd=fwdnav(18:end);

Pfwd=mat2vec_v000(fwdcov);
X=inv(eye(21)+Pfwd*PbckInv);    %%In a real implementation you must replace inv operation with cholesky based inverse. DO NOT EVER compute inverse with "inv"
W=Pfwd*X';
Y=eye(21)-W*PbckInv; %See Maybeck v2-p10

%Compute the difference between forward solution and nominal backward
%solution (i.e. obs=yb-yf)
dxDif=zeros(size(Pfwd,1),1);

%pos difference
Rappx=6378137/(sqrt(1.0-0.00669437999014*Cen_fwd(3,3)^2))+h_fwd;
vr_a=Rappx*Cen_bck(1:2,:)*Cen_fwd(3,:)';
dxDif(1:3)=[vr_a;h_fwd-h_bck];

%Vel dif
dxDif(4:6)=Vn_bck-Vn_fwd;

%Att dif
qnb_fwd=[qbn_fwd(1);-qbn_fwd(2:4)];
vr_a=quatmult_v000(qbn_bck, qnb_fwd);
mx_a=quat2dcm_v000(vr_a);
dxDif(7)=mx_a(2,3);
dxDif(8)=mx_a(3,1);
dxDif(9)=mx_a(1,2);

%imu errors
dxDif(10:end)=imu_fwd;

%Smoothed covariance
%P_smt=inv(PfwdInv+PbckInv);
P_smt=Y*Pfwd*Y'+W*PbckInv*W';

%Smoothed error estimate on the backward results
dxSmt=X*dxDif+P_smt*dzBck;

%Smoothed results
[qbn_smt, Vn_smt, Cen_smt, h_smt]=correctnav_Cen_v000(qbn_bck, Vn_bck, Cen_bck, h_bck, dxSmt(1:9), 2, 2);
imu_smt=dxSmt(10:end);

