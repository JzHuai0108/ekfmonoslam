function N = R1(rad)
% Suppose right handed coordinate frames A and B have the same x axis, 
% and rotate the A frame along the x axis by rad obtains the B frame, then
% a point in B, p_B, can be transformed from its coordinates in A frame by
% p_B = R1(rad) * p_A.
N = [1 0 0;0 cos(rad) sin(rad);0 -sin(rad) cos(rad)];