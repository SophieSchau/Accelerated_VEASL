% Returns the golden ratios for 3D imaging
function Phis = GoldenRatios3D

% Define the Fibonacci matrix as per Chan et al MRM 2009
M = [0 1 0; ...
     0 0 1; ...
     1 0 1];
 
% Derive the real eigenvector
[V,~] = eig(M);

% Normalisation is different in Chan paper, so force last component to be 1
Phis = V(1:2,1)./V(3,1);