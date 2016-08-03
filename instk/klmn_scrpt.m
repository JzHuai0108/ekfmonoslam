innoCovInv=inv(H*P*H'+R);
K=P*H'*innoCovInv;
dx=K*inno;
P=(eye(21)-K*H)*P*(eye(21)-K*H)'+K*R*K';

%Correct current states
Llh=Llh-dx(1:3);
Vn=Vn-dx(4:6);
qet=rvec2quat_v000(dx(7:9));
qbn=quatmult_v000(qet, qbn);
imuerrors=imuerrors+dx(10:end);

%smoother stuff
smo_dx_inp=smo_dx_inp+STM_upd'*H'*innoCovInv*inno(:);
smo_dP_inp=smo_dP_inp+STM_upd'*H'*innoCovInv*H*STM_upd;
STM_upd=(eye(21)-K*H)*STM_upd;