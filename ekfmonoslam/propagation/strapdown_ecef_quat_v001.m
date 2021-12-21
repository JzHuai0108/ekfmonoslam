function [qbe_new, Ve_new, ecef_new]=strapdown_ecef_quat(qbe, Ve, ecef, a, w, dt)
%Gravity (most time consuming part of ecef implementations)
Llh=ecef2geo_v000(ecef,0);
[Rn, Re, g, sL, cL, WIE_E]=geoparam_v000(Llh);
Cne=pos2Cne_v000(Llh(1), Llh(2));
ge=Cne(:,3)*g-[ecef(1:2);0]*WIE_E^2;

%Update attitude
qb=rvec2quat_v000(w*dt);
qe=rvec2quat_v000(-[0;0;WIE_E]*dt);

vr_a=quatmult_v000(qbe,qb);
qbe_new=quatmult_v000(qe,vr_a);

%%Update Velocity
%Update Vel
vel_inc1=quatrot_v000(qbe,a*dt,0);
vel_inc2=(ge+2*cross(Ve,[0;0;WIE_E])+[ecef(1:2);0]*WIE_E^2)*dt;
Ve_new=Ve+vel_inc1+vel_inc2;

%Update_pos
ecef_new=ecef+Ve*dt;

