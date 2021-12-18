function [Cbe_new, Ve_new, ecef_new]=strapdown_ecef_dcm(Cbe, Ve, ecef, a, w, dt)
%Gravity (most time consuming part of ecef implementations)
Llh=ecef2geo_v000(ecef,0);
[Rn, Re, g, sL, cL, WIE_E]=geoparam_v000(Llh);
Cne=pos2Cne_v000(Llh(1), Llh(2));
ge=-Cne(:,3)*g;

%Update attitude
rot=w*dt;
rot_norm=norm(rot);
sr_a=1-(rot_norm^2/6)+(rot_norm^4/120);
sr_b=(1/2)-(rot_norm^2/24)+(rot_norm^4/720);
mx_a=eye(3)+sr_a*skew(rot)+sr_b*skew(rot)*skew(rot);

rot=-([0;0;WIE_E])*dt;
rot_norm=norm(rot);
sr_a=1-(rot_norm^2/6)+(rot_norm^4/120);
sr_b=(1/2)-(rot_norm^2/24)+(rot_norm^4/720);
mx_b=eye(3)+sr_a*skew(rot)+sr_b*skew(rot)*skew(rot);

Cbe_new=mx_b*Cbe*mx_a;

%%Update Velocity
%Update Vel
vel_inc1=Cbe*a*dt;
vel_inc2=(-ge+2*cross(Ve,[0;0;WIE_E]))*dt;
Ve_new=Ve+vel_inc1+vel_inc2;

%Update_pos
ecef_new=ecef+Ve*dt;

