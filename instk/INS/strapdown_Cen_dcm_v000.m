%typ==0 -> nav=geodetic_frame, typ==1 -> nav=wander_frame
%although this implementation can be used as wander frame mechanization,
%this is still singular as it explicitly computes the wander angle for
%transport rate computation.
%For non-singular implementation see strapdown_wander_quat (which uses
%geowander to compute earth curvature)


function [Cbn_new, Vn_new, Cen_new, h_new]=strapdown_Cen_dcm(Cbn, Vn, Cen, h, a, w, dt, typ)

%%Geo Params
[ll, wander]=dcm2llh_v000(Cen);
Cgn=[wander(2) wander(1) 0;-wander(1) wander(2) 0;0 0 1];
llh=[ll;h];
[Rn, Re, g, sL, cL, WIE_E]=geoparam_v000(llh);

wie_n=Cen*[0; 0; WIE_E];
wen_n=Cgn*[0 1/(Re+llh(3)) 0; -1/(Rn+llh(3)) 0 0;0 -sL/cL/(Re+llh(3)) 0]*Cgn'*Vn;

if (typ==1) %wander transport rate
    wen_n(3)=0;
end

%%%Update attitude
%Part I:Body frame update
rot=w*dt;
rot_norm=norm(rot);
sr_a=1-(rot_norm^2/6)+(rot_norm^4/120);
sr_b=(1/2)-(rot_norm^2/24)+(rot_norm^4/720);
mx_a=eye(3)+sr_a*skew(rot)+sr_b*skew(rot)*skew(rot);

%Part II:Nav Frame update
rot=-(wen_n+wie_n)*dt;
rot_norm=norm(rot);
sr_a=1-(rot_norm^2/6)+(rot_norm^4/120);
sr_b=(1/2)-(rot_norm^2/24)+(rot_norm^4/720);
mx_b=eye(3)+sr_a*skew(rot)+sr_b*skew(rot)*skew(rot);

Cbn_new=mx_b*Cbn*mx_a;

%%%Update Velocity
vel_inc1=(Cbn*(a*dt))+[0;0;g]*dt;
vel_inc2=(cross(Vn,2*wie_n+wen_n))*dt;
Vn_new=Vn+vel_inc1+vel_inc2;

%%%Update Cen (position + wander_angle)
rot=-wen_n*dt;
rot_norm=norm(rot);
sr_a=1-(rot_norm^2/6)+(rot_norm^4/120);
sr_b=(1/2)-(rot_norm^2/24)+(rot_norm^4/720);
mx_b=eye(3)+sr_a*skew(rot)+sr_b*skew(rot)*skew(rot);
Cen_new=mx_b*Cen;

%update height
h_new=h-Vn(3)*dt;
