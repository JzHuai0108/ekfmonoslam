%wander=[sin(alpha) cos(alpha)]
%R=[Rn;Re];
%If you use llh parametrization, you have to update wander angle manually.
%Note than, when you use Cen as position parametrization, wander angle is
%automatically updated by setting the z-component of transport rate to "0"

%I had used this in a large heading error implementation assuming that the
%wander represents the unknown large heading.

%Under small angle assuption, do not use this implementation.

%Note:This implementation calls geoparam_v000. Therefore, this is still
%singular at the poles. See strapdown_wander_quat_v000 for a non-singular
%implementation (which is essentially the same as this except it calls
%geoparam_v001 and uses "Cen" for position parametrization.) (Obviously, it is not
%possible to obtain a singular implementation with Llh parametrization: Latitude
%is not defined at the poles).

%(For an example usage see: example_largeheading)

function [Cbn_new, Vn_new, llh_new, wander_new]=strapdown_wander_dcm(Cbn, Vn, llh, wander, a, w, dt)

wander=wander(:); %just to make it a column matrix;

%Geo Parameters
[Rn, Re, g, sL, cL, WIE_E]=geoparam_v000(llh);
Cgn=[wander(2) wander(1) 0;-wander(1) wander(2) 0;0 0 1];

wen_n=Cgn*[0 1/(Re+llh(3)) 0; -1/(Rn+llh(3)) 0 0;0 0 0]*Cgn'*Vn;
wie_n=Cgn*[WIE_E*cL; 0; -WIE_E*sL];

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

%%%Update position
llh_new=llh+[1/(Rn+llh(3)) 0 0;0 1/cL/(Re+llh(3)) 0;0 0 -1]*Cgn'*Vn*dt;

%%%update wander angle
alpha_dot=Vn(2)/(Re+llh(3))*(sL/cL);
wander_new=wander+[wander(2);-wander(1)]*alpha_dot*dt;