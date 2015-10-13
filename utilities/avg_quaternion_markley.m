
% Averaging Quaternions
% Since quaternions are not regular vectors, but rather representations 
% of orientation, an average quaternion cannot just be obtained by taking 
% a weighted mean. This function implements the work done by paper by 
% F. Landis Merkley to calculate the average quaternion. The algorithm 
% explained by F. Landis Markley at: 
% http://www.acsu.buffalo.edu/~johnc/ave_quat07.pdf
% For this particular implementation, I would also like to reference Mandar
% Harshe:
% http://www-sop.inria.fr/members/Mandar.Harshe/knee-joint/html/index.html
%
% This algorithm is compared by rotqrmean from VoiceBox and found to 
% produce quite similar results, yet it is more elegant, much simpler to 
% implement and follow. (Though, there might be difference in signs)
%
% Usage : 
% Q is an Mx4 matrix, where each row stores a quaternion to be averaged.
% In return, the function outputs Qavg, which is a single quaternion
% corresponding to the average.
%
% Tolga Birdal

function [Qavg]=avg_quaternion_markley(Q)

% Form the symmetric accumulator matrix
A=zeros(4,4);
M=size(Q,1);

for i=1:M
    q=Q(i,:)';
    A=q*q'+A; % rank 1 update
end

% scale
A=(1.0/M)*A;

% Get the eigenvector corresponding to largest eigen value
[Qavg, Eval]=eigs(A,1);

end
