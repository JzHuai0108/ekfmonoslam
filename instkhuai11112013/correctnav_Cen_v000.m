%Corrects the navigation states using the Kalman generated error values
%if size(att)=[4,1], it assumes att=qbn else it assumes att=Cbn
%vel=Vn
%dx=corrections where position errors are defined as dr^n (order:[pos, vel, att])
%modtype: 1=PHI model, 2=PSI model
%mechtype: 1:wander-azimuth 2:local-geodetic

function [att_new, Vn_new, Cen_new, h_new]=correctnav_Cen(att, Vn, Cen, h, dx, modtype, mechtype)

%Correct the position
Fc=geoparam_v001(mechtype, Cen(:,3), h, []);
theta=Fc*[dx(2);-dx(1);0];
Ccn=rot2dcm_v000(theta);
Cen_new=Ccn*Cen;
h_new=h+dx(3);  %note that dr_z and dh_z are in opposite direction. That is why we add the error.

%velocity correction
Vn_new=Vn-dx(4:6);

%attitude correction
if (size(att,1)==3) %attitude is represented as DCM
    Cet=rot2dcm_v000(dx(7:9));   %Erroneous to true navigation frame:DCM=(I+S(ang))
    att_new=Cet*att;
else %attitude is quaternion
    %Method I
    qet=rvec2quat_v000(dx(7:9));
    att_new=quatmult_v000(qet, att);
    
%     %Method II:
%     %Method I involves some trigonometric functions. Here is a
%     %algebraic correction method (See: Chung 1996-Eq:16) (Special thanks to Kaygisiz)
%     R=[-q(2) -q(3) -q(4);q(1) q(4) -q(3);q(-4) q(1) q(2);q(3) q(-2) q(1)];
%     dq=-0.5*R*dx(7:9);
%     att_new=att_new-dq;
%     att_new=att_new/norm(att_new);  %optional
end


if (modtype==2) %dx represents the psi model errors.
    %In this case, when we change Cen, we implicitly change the "c" frame.
    %Therefore, we also need to rotate attitude and velocity accordingly for the new (updated) "c" frame.
    
    Vn_new=Ccn*Vn_new;
    
    if (size(att,1)==3)
        att_new=Ccn*att_new;
    else
        att_new=quatmult_v000(rvec2quat_v000(theta), att_new);
    end
end