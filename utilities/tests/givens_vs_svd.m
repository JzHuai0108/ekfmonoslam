% compare 3 methods to compute null space basis vectors: svd, 
% matlab null which is said to use svd, and qr decomposition
% qr decomposition may be done with givens or householder rotations

trials=1e1;
timers=zeros(3,1);
for shot=1:trials
Hf=rand(10,3);

tic
[u,d,v]= svd(Hf');
n1= v(:, 4:end);
timers(1) = timers(1) + toc;
% assert(length(find(d>0))==3);
%svd to get nullspace basis, https://inst.eecs.berkeley.edu/~ee127a/book/login/l_svd_mat_prop.html
disp(n1);

tic
n2=null(Hf');
timers(2) = timers(2) + toc;
disp(n2);

tic
% compute nullspace with qr decomposition
% http://stackoverflow.com/questions/2181418/computing-the-null-space-of-a-matrix-as-fast-as-possible
[q, r]= qr(Hf);
n3=q(:,4:end);
timers(3) = timers(3) + toc;
disp(n3);
end
disp('times consumed for computing nullspace by svd, matlab null, and qr');
disp(timers)


% To project a matrix or vector onto the nullspace of another matrix is
% useful in eliminating unnecessary variables. e.g., r= H_x*X+H_f*p +n
% where H_f is of size 2M*3. p is of size 3*1. To remove p, we need to
% project r and H_x onto the nullspace of H_f', Q_2.
% Note H_f =[Q_1, Q_2][R; 0]= Q[R; 0]
% 
% To accomplish this, we can successively apply Givens rotations on 
% [r, H_x, H_f] until it becomes [Q'*r, Q'*H_x, [R; 0]]. And remove its first
% 3 rows, we get [Q_2'*r, Q_2'*H_x, 0].
% 
% This trick is used in Mourikis, Anastasios, and Stergios Roumeliotis. 
% "A multi-state constraint Kalman filter for vision-aided inertial navigation." 
% Robotics and Automation, 2007 IEEE International Conference on. IEEE, 2007.
% How to apply Givens rotations can be found in 
% G. Golub and C. van Loan, Matrix computations. The Johns Hopkins
% University Press, London, 1996.

trials=1e2;
timers=zeros(2,1);
% For small matrices, svd outperforms qr +svd; for large thin tall
% matrices, svd is a little less efficient than qr +svd
m =10; %3000;
n =5; %1200;
for shot=1:trials
assert(m>n);
A= rand(m, n);
tic
[~, ~, v] = svd(A);
x1= v(:, end);
timers(1) = timers(1) + toc;

tic
[~, r] = qr(A);
Th= r(1:n, 1:n);
[~, ~, v2]= svd(Th);
x2= v2(:, end);
timers(2) = timers(2) + toc;

end

disp('times consumed for solving least squares by svd, and qr+svd');
disp(timers)

