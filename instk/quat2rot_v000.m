%Convert quaternion into rotation vector
%assumes that sin(theta/2) is +ve
function rot=quat2rot(quat)
sr_a=sqrt(quat(2)^2+quat(3)^2+quat(4)^2);
if (sr_a==0)
    rot=zeros(3,1);
else
    theta=2*atan2(sr_a, quat(1));
    rot=theta*[quat(2);quat(3);quat(4)]/sin(theta/2);
end
