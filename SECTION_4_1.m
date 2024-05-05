% This scripts tests the correctness of the numerically computed 
% VEM local matrices for the unit square

% Define array P containing coordinates of vertexes of F
P = [0 0 0; 0 1 0; 1 1 0; 1 0 0];
% F is star-shaped w.r.t. the following point P0
P0 = [.5 .5 0];
% Create the unit square F and compute its local VEM matrices
F = getLocalMatrices(element2d(P, false, false, [], P0));
% Closed-form local VEM matrices
K = (4*eye(4)-ones(4))/4;
M = [17 -9 13 -9;
     -9 17 -9 13;
     13 -9 17 -9;
     -9 13 -9 17]/48;
C = [5 3 1 3;
     3 5 3 1;
     1 3 5 3;
     3 1 3 5]/48;
% Evaluate maximum error of numerically computed local matrices
norm(K(:) - F.K(:),'inf')
norm(C(:) - F.C(:),'inf')
norm(M(:) - F.M(:),'inf')