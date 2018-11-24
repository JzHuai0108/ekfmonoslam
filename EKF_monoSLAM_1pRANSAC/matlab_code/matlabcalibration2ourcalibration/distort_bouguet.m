function distorted = distort_bouguet(undistorted,kc,fc,cc)

xu = ( undistorted(1,:) - cc(1) )/fc(1);
yu = ( undistorted(2,:) - cc(2) )/fc(2);

    k1 = kc(1);
    k2 = kc(2);
    k3 = kc(5);
    p1 = kc(3);
    p2 = kc(4);
    
r_2 = sum([xu;yu].^2);
k_radial =  1 + k1 * r_2 + k2 * r_2.^2 + k3 * r_2.^3;

dx = [2*p1*xu.*yu + p2*(r_2+2*xu.^2); 2*p2*xu.*yu + p1*(r_2+2*yu.^2)];

distorted = [xu;yu].*[1 + k1 * r_2 + k2 * r_2.^2 + k3 * r_2.^3; 1 + k1 * r_2 + k2 * r_2.^2 + k3 * r_2.^3] + dx;

distorted = [ distorted(1,:)*fc(1) + cc(1);  distorted(2,:)*fc(2) + cc(2)  ];