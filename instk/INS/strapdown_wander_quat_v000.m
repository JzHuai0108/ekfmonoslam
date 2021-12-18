%Wander-azimuth strapdown mechanization with quaternion attitude
%implementation. Geodetic position is updated using Cen.

%In contrast to strapdown_wander_dcm, this is a non-singular implementation.
%Yes, this will work in poles without any problems.

%For an example see: example_wander.m

function [qbn_new, Vn_new, Cen_new, h_new]=strapdown_wander_quat(qbn, Vn, Cen, h, velinc, anginc, dt)
%% Compute the transport rate and geodetic parameters (pz=0)
[Fc, wen_n, wie_n, g]=geoparam_v001(1, Cen(:,3), h, Vn); %non-singular implementation

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
%Note that we have here both Vn and Vn_new. Why not use both to update position?
%Savage loves combining old and new parameters as a trapezoidal integral.
%Let us follow his advice here. In fact, I could have done this for each
%step above.
Vmid=(Vn+Vn_new)/2;
wen_n=Fc*[Vmid(2);-Vmid(1);0];
Cn=rot2dcm_v000(-wen_n*dt);
Cen_new=Cn*Cen;

h_new=h-Vmid(3)*dt;