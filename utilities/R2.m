function N = R2(rad)
% Suppose right handed coordinate frames A and B have the same y axis, 
% and rotate the A frame along the y axis by rad obtains the B frame, then
% a point in B, p_B, can be transformed from its coordinates in A frame by
% p_B = R2(rad) * p_A.
N = [cos(rad) 0 -sin(rad);0 1 0;sin(rad) 0 cos(rad)];