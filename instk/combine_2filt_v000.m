%Combine forward and backward solutions

%The input arguments is highly dependent on the implementation. In the future I
%will probably change script with some other implementation independent
%version.

%dzBck=PbckInv*dx_bckwrd
function [P_smt, Cen_smt, h_smt, Vn_smt, qbn_smt, imu_smt, dxDif]=combine_2filt(fwdnav, fwdcov, Cen_bck, h_bck, Vn_bck, qbn_bck, dzBck, PbckInv)

%Forward filter results
Cen_fwd=reshape(fwdnav(1:9),[3,3]);
h_fwd=fwdnav(10); %elipsoid height
Vn_fwd=fwdnav(11:13);
qbn_fwd=fwdnav(14:17);
imu_fwd=fwdnav(18:end);

Pfwd=mat2vec_v000(fwdcov);
%%IMPORTATNT::In a real implementation you must replace inv operation with cholesky based inverse. DO NOT EVER compute inverse with "inv"
X=inv(eye(size(Pfwd))+Pfwd*PbckInv); 
Wt=X*Pfwd;      %Note that (Maybeck, v2p10) uses the transpose of W. I did not understand why he prefers the transpose. (It seems it complicates the results. Could be because of some numeric properties?)

%Compute the difference between the forward solution and the nominal backward
%solution (we are going to estimate the errors on nominal forward solution
%using backward. We choose this approach because PfwdInv*dxFwd=0 as dxFwd=0)
dxDif=zeros(size(Pfwd,1),1);

%pos difference
Rappx=6378137/(sqrt(1.0-0.00669437999014*Cen_fwd(3,3)^2))+h_fwd;
vr_a=Rappx*Cen_fwd(1:2,:)*Cen_bck(3,:)';
dxDif(1:3)=[vr_a;h_bck-h_fwd];

%Vel dif
dxDif(4:6)=Vn_fwd-Vn_bck;

%Att dif
qnb_bck=[qbn_bck(1);-qbn_bck(2:4)];
vr_a=quatmult_v000(qbn_fwd, qnb_bck);
mx_a=quat2dcm_v000(vr_a);
dxDif(7)=mx_a(2,3);
dxDif(8)=mx_a(3,1);
dxDif(9)=mx_a(1,2);

%imu errors
dxDif(10:end)=imu_fwd;

%Smoothed covariance
%P_smt=inv(PfwdInv+PbckInv);
P_smt=X*Pfwd*X'+Wt*PbckInv*Wt';

%Smoothed error estimate on the forward results
dxSmt=P_smt*(PbckInv*dxDif+dzBck);

%Smoothed results
[qbn_smt, Vn_smt, Cen_smt, h_smt]=correctnav_Cen_v000(qbn_fwd, Vn_fwd, Cen_fwd, h_fwd, dxSmt(1:9), 2, 2);
imu_smt=imu_fwd-dxSmt(10:end);

