%Quaternion based strapdown implementation on geodetic frame.
%Position is represented as lat, long and height

function [qbn_new, Vn_new, llh_new]=strapdown_ned_quat(qbn, Vn, llh, vel_inc, ang_inc, dt)
%Geo Parameters
[Rn, Re, g, sL, cL, WIE_E]=geoparam_v000(llh);

wen_n=[Vn(2)/(Re+llh(3)); -Vn(1)/(Rn+llh(3)); -Vn(2)*sL/cL/(Re+llh(3))];
wie_n=[WIE_E*cL; 0; -WIE_E*sL];

%% Update the velocity
vel_inc1=quatrot_v000(qbn,vel_inc,0);
vel_inc2=(cross(Vn,2*wie_n+wen_n)+[0;0;g])*dt;
Vn_new=Vn+vel_inc1+vel_inc2;

%Now I can use the middle velocity value (This does not make any
%difference, but ,just for the sake of different flavour, I am going to use it)
Vn_mid=(Vn+Vn_new)/2;

%%update the attitude
%Body frame updates
qb=rvec2quat_v000(ang_inc);
qbn_new=quatmult_v000(qbn,qb);

%Navigation frame updates
wen_n=[Vn_mid(2)/(Re+llh(3)); -Vn_mid(1)/(Rn+llh(3)); -Vn_mid(2)*sL/cL/(Re+llh(3))];    %just another useless modification.
qn=rvec2quat_v000(-(wen_n+wie_n)*dt);
qbn_new=quatmult_v000(qn,qbn_new);

%%Update position
llh_new=llh+[1/(Rn+llh(3)) 0 0;0 1/cL/(Re+llh(3)) 0;0 0 -1]*Vn_mid*dt;
