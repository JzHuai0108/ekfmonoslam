function R = eul2rotm_oxford(r, p, yaw)
% copied from matlab/SE3MatrixFromComponents.m oxford robotcar dataset.
  R_x = [ 
    1, 0, 0;
    0, cos(r), -sin(r);
    0, sin(r), cos(r) ];
  
  R_y = [ 
    cos(p), 0, sin(p);
    0, 1, 0;
    -sin(p), 0, cos(p) ];
  
  R_z = [
    cos(yaw), -sin(yaw), 0;
    sin(yaw), cos(yaw), 0;
    0, 0, 1];
  
  R = R_z * R_y * R_x;