%nav formulas on a (non rotating flat earth)
%vtyp==0->n frame vel, vtyp==1->b frame vel
function [Cbn_new, Vx_new, pos_new]=strapdown_pln_dcm(Cbn, Vx, pos, a, w, g, dt, vtyp)

%%%Update attitude
%Part I:Body frame update
rot=w*dt;
rot_norm=norm(rot);
sr_a=1-(rot_norm^2/6)+(rot_norm^4/120);
sr_b=(1/2)-(rot_norm^2/24)+(rot_norm^4/720);
mx_a=eye(3)+sr_a*skew(rot)+sr_b*skew(rot)*skew(rot);

Cbn_new=Cbn*mx_a;

%%Update Velocity and position
if (vtyp==0) %vel in n frame
    vel_inc=(Cbn*(a*dt))+[0;0;g]*dt;
    Vx_new=Vx+vel_inc;
    
    pos_new=pos+Vx*dt;
elseif (vtyp==1) %vel in b frame
    vel_inc1=(a+(Cbn'*[0;0;g]))*dt;
    vel_inc2=(cross(Vx,w))*dt;
    Vx_new=Vx+vel_inc1+vel_inc2;
   
    pos_new=pos+Cbn*Vx*dt;
end