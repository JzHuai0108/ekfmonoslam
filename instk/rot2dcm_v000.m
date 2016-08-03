%"a" is rotated with "rot" to obtain "b".
%dcm=Cba (from b to a)

function dcm=rot2dcm(rot)
rot_norm=norm(rot);
if (rot_norm>0)
    sr_a=sin(rot_norm)/rot_norm;
    sr_b=(1-cos(rot_norm))/rot_norm^2;
    
    dcm=eye(3)+sr_a*skew(rot)+sr_b*skew(rot)*skew(rot);
else
    dcm=eye(3);
end


