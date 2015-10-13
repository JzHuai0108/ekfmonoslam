function [Cbn_new, Vb_new, llh_new]=strapdown_bd_dcm(Cbn, Vb, llh, a, w, dt)
%Geo Parameters
[Rn, Re, g, sL, cL, WIE_E]=geoparam_v000(llh);

%N-frame Vel
Vn=Cbn*Vb;

%N-frame rotation
wen_n=[Vn(2)/(Re+llh(3)); -Vn(1)/(Rn+llh(3)); -Vn(2)*sL/cL/(Re+llh(3))];
wie_n=[WIE_E*cL; 0; -WIE_E*sL];

%Update attitude
rot=w*dt;
rot_norm=norm(rot);
sr_a=1-(rot_norm^2/6)+(rot_norm^4/120);
sr_b=(1/2)-(rot_norm^2/24)+(rot_norm^4/720);
mx_a=eye(3)+sr_a*skew(rot)+sr_b*skew(rot)*skew(rot);

rot=-(wen_n+wie_n)*dt;
rot_norm=norm(rot);
sr_a=1-(rot_norm^2/6)+(rot_norm^4/120);
sr_b=(1/2)-(rot_norm^2/24)+(rot_norm^4/720);
mx_b=eye(3)+sr_a*skew(rot)+sr_b*skew(rot)*skew(rot);

Cbn_new=mx_b*Cbn*mx_a;

%%Update Velocity
vel_inc1=(a+(Cbn'*[0;0;g]))*dt;
vel_inc2=(cross(Vb,(Cbn'*wie_n)+w))*dt;
Vb_new=Vb+vel_inc1+vel_inc2;

%Update_pos
llh_new=llh+[1/(Rn+llh(3)) 0 0;0 1/cL/(Re+llh(3)) 0;0 0 -1]*Vn*dt;
