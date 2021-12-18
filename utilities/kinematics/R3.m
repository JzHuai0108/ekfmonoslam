function N = R3(rad)
% Suppose right handed coordinate frames A and B have the same z axis, 
% and rotate the A frame along the z axis by rad obtains the B frame, then
% a point in B, p_B, can be transformed from its coordinates in A frame by
% p_B = R3(rad) * p_A.
N = [cos(rad) sin(rad) 0;-sin(rad) cos(rad) 0;0 0 1];