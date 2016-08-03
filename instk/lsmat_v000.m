function [LSM LSM_R] =lsmat(H,R)
R_inv=inv(R);
LSM_R=inv(H'*R_inv*H);
LSM=LSM_R*H'*R_inv;