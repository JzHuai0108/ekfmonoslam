function xyd_xyn=xydistort_by_xynormalized(xyn, cam)
k1=cam.kc(1);
k2=cam.kc(2);
p1=cam.kc(3);
p2=cam.kc(4);
k3=cam.kc(5);
xn=xyn(1);
yn=xyn(2);
xyd_xyn=zeros(2);
% xyd to xyn
xyd_xyn(1,1)=k2*(xn^2 + yn^2)^2 + k3*(xn^2 + yn^2)^3 + 6*p2*xn + 2*p1*yn + xn*(2*k1*xn ...
+ 4*k2*xn*(xn^2 + yn^2) + 6*k3*xn*(xn^2 + yn^2)^2) + k1*(xn^2 + yn^2) + 1;
% xd to yn
xyd_xyn(1,2)=2*p1*xn + 2*p2*yn + xn*(2*k1*yn + 4*k2*yn*(xn^2 + yn^2) + 6*k3*yn*(xn^2 + yn^2)^2); 
% yd to (xn, yn) 
xyd_xyn(2,1)=2*p1*xn + 2*p2*yn + yn*(2*k1*xn + 4*k2*xn*(xn^2 + yn^2) + 6*k3*xn*(xn^2 + yn^2)^2);
xyd_xyn(2,2)=k2*(xn^2 + yn^2)^2 + k3*(xn^2 + yn^2)^3 + 2*p2*xn + 6*p1*yn + yn*(2*k1*yn ...
+ 4*k2*yn*(xn^2 + yn^2) + 6*k3*yn*(xn^2 + yn^2)^2) + k1*(xn^2 + yn^2) + 1;
 