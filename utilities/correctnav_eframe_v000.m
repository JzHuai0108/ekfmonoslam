%Corrects the navigation states using the Kalman generated error values
%if size(att)=[4,1], it assumes att=qbn else it assumes att=Cbn
%vel=Vn
%dx=corrections where position errors are defined as dr^n (order:[pos, vel, att])
%delta_x(k)=K(delta_y(k)-h(delta_x(k|k-1)))=K*delta_y(k)

function [att_new, Ve_new, ecef_new]=correctnav_eframe_v000(att, Ve, ecef, dx)

%Correct the position
ecef_new=ecef-dx(1:3);

%velocity correction
Ve_new=Ve-dx(4:6);

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
