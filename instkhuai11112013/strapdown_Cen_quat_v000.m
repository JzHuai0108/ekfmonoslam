%Navigation frame=Local geodetic Frame.
%Position is mechanized in Cen (essentially Ceg) and h.
%Attitude is mechanized in quat

%Compare this with strapdown_Cen_dcm and strandown_wander_quat. Both of
%these strapdowns implement wander frame mechanization. However, only the
%2nd one is non-singular. Furthermore, also note that the only difference
%between this script and the strandown_wander_quat is the 1st argument of
%geoparam function.

function [qbn_new, Vn_new, Cen_new, h_new]=strapdown_Cen_quat(qbn, Vn, Cen, h, velinc, anginc, dt)
%Compute the transport rate of geodetic frame (also the gravity)
[Fc, wen_n, wie_n, g]=geoparam_v001(2, Cen(:,3), h, Vn); %Compute the curvature matrix for local geodetic frame

%% Update the velocity
vel_inc1=quatrot_v000(qbn,velinc,0); %Note: 'a' is assumed to be output of a sculling module.
vel_inc2=(cross(Vn,2*wie_n+wen_n)+[0;0;g])*dt;
Vn_new=Vn+vel_inc1+vel_inc2;

%%update the attitude
%Body frame updates
qb=rvec2quat_v000(anginc);    %Note: 'w' is assumed to be the output of coning module
qbn_new=quatmult_v000(qbn,qb);

%Navigation frame updates
qn=rvec2quat_v000(-(wen_n+wie_n)*dt);
qbn_new=quatmult_v000(qn,qbn_new);

%%Position update
Vmid=(Vn+Vn_new)/2;
wen_n=Fc*[Vmid(2);-Vmid(1);0];
Cn=rot2dcm_v000(-wen_n*dt);
Cen_new=Cn*Cen;

h_new=h-Vmid(3)*dt;