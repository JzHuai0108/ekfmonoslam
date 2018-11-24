function [xd, dxd_dx] = apply_fisheye_distortion(x,k)

%apply_fisheye_distortion.m
%
%[x] =  apply_fisheye_distortion(xd,k)
%
%Apply the fisheye distortions
%
%INPUT: x: undistorted (normalized) point coordinates in the image plane (2xN matrix)
%       k: Fisheye distortion coefficients (5x1 vector)
%
%OUTPUT: xd: distorted (normalized) point coordinates in the image plane (2xN matrix)
%        dxd_dx : jacobian (2x2xN matrix)

r = sqrt(x(1,:).^2 + x(2,:).^2);

theta = atan(r);
theta_d = theta .* (1 + k(1)*theta.^2 + k(2)*theta.^4 + k(3)*theta.^6 + k(4)*theta.^8);

scaling = ones(1,length(r));

ind_good = find(r > 1e-8);

scaling(ind_good) = theta_d(ind_good) ./ r(ind_good);

xd = x .* (ones(2,1)*scaling);

if (nargout > 1)
  n = size(x,2);
  dr_dx = x;  % 2xn
  dr_dx(:,ind_good) = dr_dx(:,ind_good) ./ repmat(r(ind_good), [2,1]);

  theta_d_r2 = scaling;
  theta_d_r2(ind_good) = theta_d_r2(ind_good) ./ r(ind_good);

  dtheta_dr = 1 ./ (1 + r.^2);
  dtheta_d_dtheta = 1 + 3*k(1)*theta.^2 + 5*k(2)*theta.^4 + ...
      7*k(3)*theta.^6 + 9*k(4)*theta.^8;
  dtheta_d_dr_r = dtheta_d_dtheta .* dtheta_dr;  % 1xn
  dtheta_d_dr_r(ind_good) = dtheta_d_dr_r(ind_good) ./ r(ind_good);

  dxd_dx = zeros(2,2,n);
  for i = 1:n
    dxd_dx(:,:,i) = scaling(i)*eye(2,2) + x(:,i)*(dtheta_d_dr_r(i) - ...
                    theta_d_r2(i))*dr_dx(:,i)';
  end
end
