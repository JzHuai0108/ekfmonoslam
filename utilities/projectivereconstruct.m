function [Int, Ext]=projectivereconstruct(Ab)
% recover Intrinsic and Extrinsic Parameters from the projection matrix
% also see computer vision a modern approach section 3.2 first edition by
% Forsyth and Ponce
% written by Jianzhu Huai, March 2014
% Ab is the Projection Matrix, 3 X 4 matrix
% example
% K=[1000, 0, 501; 0, 900, 460; 0,0, 1];
% M=[R3(-.2)*R2(-0.5)*R1(0.9), rand(3,1)]
% Ab=K*M;
% [Int, Ext]=projectivereconstruct(Ab)
a1 = Ab(1,1:3); a2 = Ab(2,1:3); a3 = Ab(3,1:3);b = Ab(:,4);
% Computations of Parameters
sqa3=(a3*a3');
ep = -1; % ep = +1 or -1
r3 = ep*a3/sqrt(sqa3);
u0 = a1*a3'/sqa3;
v0 = a2*a3'/sqa3;
theta = acos(-(dot(cross(a1,a3),cross(a2,a3)))/(norm(cross(a1,a3),2)*norm(cross(a2,a3),2)));
alpha = norm(cross(a1,a3),2)*sin(theta)/sqa3;
beta = norm(cross(a2,a3),2)*sin(theta)/sqa3;
r1 = (1/norm(cross(a2,a3),2))*cross(a2,a3);
r2 = cross(r3,r1);
K = [ alpha -alpha*cot(theta) u0; 0 beta/sin(theta) v0; 0 0 1];
t =ep*inv(K)*b/sqrt(sqa3);
if t(3)<0 
    t=-t;
    r3=-r3;
    r2=-r2;
end
R = [r1;r2;r3];
% Final Results
Ext = [R t]; % Extrinsic Parameters
Int = K; % Intrinsic Parameters
M = Int*Ext; % Projection Matrix
end