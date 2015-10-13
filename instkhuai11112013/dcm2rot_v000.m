%dcm (Cba) to rotation vector (r^a)
%Cba: B to A transformation matrix
%r^a=Rotation vector defined in A. (Rotate A around r with mag(r) to obtain B)
%the magnitude of the returning rotation vector will be between [0 and PI]

%This is the implementaion of savage(3-31).
function rvec=dcm2rot(dcm)

sinPHI=0.5*sqrt((dcm(3,2)-dcm(2,3))^2+(dcm(1,3)-dcm(3,1))^2+(dcm(2,1)-dcm(1,2))^2);
cosPHI=0.5*(dcm(1,1)+dcm(2,2)+dcm(3,3)-1);

PHI=atan2(sinPHI,cosPHI);

if (cosPHI>=0)
    if (sinPHI==0)
        F=1;
    else
        F=PHI/sinPHI;
    end
    rvec=0.5*F*[dcm(3,2)-dcm(2,3);dcm(1,3)-dcm(3,1);dcm(2,1)-dcm(1,2)];
else
    sr_a=1-cosPHI;
    mu1=sqrt((dcm(1,1)-1)/sr_a+1);
    mu2=sqrt((dcm(2,2)-1)/sr_a+1);
    mu3=sqrt((dcm(3,3)-1)/sr_a+1);
    if (mu1>=mu2 && mu1>=mu3)
       u1=mu1*sign(dcm(3,2)-dcm(2,3));
       u2=1/2/u1*(dcm(1,2)+dcm(2,1))/sr_a;
       u3=1/2/u1*(dcm(1,3)+dcm(3,1))/sr_a;
    elseif (mu2>=mu3)
       u2=mu2*sign(dcm(1,3)-dcm(3,1));
       u3=1/2/u2*(dcm(2,3)+dcm(3,2))/sr_a;
       u1=1/2/u2*(dcm(1,2)+dcm(2,1))/sr_a; 
    else
       u3=mu3*sign(dcm(2,1)-dcm(1,2));
       u1=1/2/u3*(dcm(1,3)+dcm(3,1))/sr_a;
       u2=1/2/u3*(dcm(2,3)+dcm(3,2))/sr_a;
    end
    rvec=PHI*[u1;u2;u3];
end