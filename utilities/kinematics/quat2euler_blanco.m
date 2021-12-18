function eul=quat2euler_blanco(q)
% implement the method described in a tutorial to SE(3) transform by Blanco
% q=[w,x,y,z]', eul=[roll, pitch, yaw]' in accord with Jekeli 2000 and Blanco
% it seems Blanco's method is a little preciser than that in Voicebox
% test cases: 
% quat=[0.5,-0.5,0.5,0.5]'; % two possible euler angle solutions
% quat=[0.5,0.5,-0.5,0.5]'; % two possible euler angle solutions
% quat= rand(4,1); quat=quat/norm(quat);
% eul1=rotqr2eu('xyz', quat); % using voicebox functions
% eul2=quat2euler_blanco(quat);
% roteu2ro('xyz', eul1)-rotqr2ro(quat)
% roteu2ro('xyz', eul2)-rotqr2ro(quat)
eul=zeros(3,1);
delta= q(1)*q(3)-q(2)*q(4);
eps=1e-8;
if(0.5-abs(delta)<eps)
    if(delta>0)
        eul(3)=-2*atan2(q(2), q(1));
        eul(2)=pi/2;
    else
        eul(3)=2*atan2(q(2), q(1));
        eul(2)=-pi/2;
    end
else
eul(3)=atan2(2*(q(1)*q(4)+q(2)*q(3)), 1-2*(q(3)^2+q(4)^2));
eul(2)=asin(2*delta);
eul(1)=atan2(2*(q(1)*q(2)+q(3)*q(4)), 1-2*(q(2)^2+q(3)^2));
end
end