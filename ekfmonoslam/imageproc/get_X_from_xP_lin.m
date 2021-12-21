function X=get_X_from_xP_lin(xn,Pcomp)
%get_X_from_xP_lin  Estimation of 3D point from image matches and camera matrices, linear.
%  return a 3D point X (column 4-vector, homogeneous, last element 1)
%   from its projections in K images xn (2-by-K matrix) and camera matrices
% P (K-cell of 3-by-4 matrices) by minimizing algebraic distance.
%   xn 3 x K matrix, homogeneous normalized image coordinates,
% last rows of 1, the other values usu. falls between -0.5 to 0.5
K=size(xn,2);
A = zeros(2*K, 4);
for k = 1:K
  A(2*k+(-1:0),:)=[xn(1,k)*Pcomp(3,:,k)-Pcomp(1,:,k); xn(2,k)*Pcomp(3,:,k)-Pcomp(2,:,k)];
end
% A = normx(A')';
[dummy,dummy,X] = svd(A,0);
X = X(:,end);
X=X/X(4);