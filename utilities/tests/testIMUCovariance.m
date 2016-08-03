function testIMUCovariance()
% see what happens to covariance after many propagations
dt=0.1;
qa=(80*1.0e-5)^2;
qg=(0.03*pi/180)^2;
gravs=9.8;
Pis=zeros(9,floor(450/dt)+1);
P0=eye(9)*1e-3;
Pis(:, 1)=sqrt(diag(P0));
P1=zeros(9);
for i=2:size(Pis,2)
    Phi=eye(9);
    Phi(1:3,4:6)=eye(3)*dt;
    covu=zeros(9);
    covu(4:6,4:6)=eye(3)*qa*dt;
    covu(7:9,7:9)=eye(3)*qg*dt;
    P1=Phi*P0*Phi'+covu;
    P0=P1;
    Pis(:, i)=sqrt(diag(P1));
end
plot((0:1:size(Pis,2)-1)*dt, Pis(1:3,:));
legend('x','y','z');
end