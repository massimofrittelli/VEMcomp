function [KB,MB,KS,MS,ES,R] = matricesTet(PB, EB)

% DESCRIPTION - Generates BSFEM matrices for a 3D tetrahedral mesh. Source:
% https://it.mathworks.com/matlabcentral/answers/678418-calculating-tetrahedral-stiffness-matrix-for-fem
% after suitable testing and adaptations.

% INPUTS:
% - PB is a NPBx3 array containing the coordinates of the nodes by rows
% - EB is a NEBx4 array containing the tetrahedral elements by rows

% OUTPUTS:
% - KB, MB: stiffness and mass matrices in the bulk
% - KS, MS: stiffness and mass matrices on the surface
% - ES: an array containing the surface elements by rows
% - R: the reduction matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Number of nodes and elements in the bulk
NPB = size(PB,1);
NEB = size(EB,1);

% Mass matrix on reference tetrahedron in 3d
MrefB = (ones(4) + eye(4))/120;

%create sparse matrix by storing (i,j) indices and the values
inds = 0;
Si = zeros(16*NEB,1); %ith indices of stiffness matrix
Sj = zeros(16*NEB,1); %jth indices of stiffness matrix
Sk = zeros(16*NEB,1); %value of (i,j) element of stiffness matrix
Mk = zeros(16*NEB,1); %value of (i,j) element of mass matrix
%4 indices vs 4 indices
ii = 1:4; ji = 1:4; [ii, ji] = ndgrid(ii,ji);
ii = ii(:); ji = ji(:);

% Loop through every element in the mesh - would be nice to do this without
% looping - or this maybe why it's best to use a compiled language
for i = 1:NEB
    
    eleNodes = EB(i,1:4);
    %uses simple integration of linear terms method
    B = [PB(eleNodes,:)'; [1, 1, 1, 1]];
    C = inv(B);
    ve = abs(det(B)); % 6 times the signed volume of the tetrahedral element
    
    Kl2 = C(:,1:3)*C(:,1:3)'*ve/6;
    
    %store all 16 terms
    inds = inds(end) + (1:16);
    Si(inds) = eleNodes(ii);
    Sj(inds) = eleNodes(ji);
    Sk(inds) = Kl2(:);
    Mk(inds) = MrefB(:)*ve;
end
%Build stiffness matrix - MATLABs sparse operator naturally adds together
%all terms with the same i and j - i.e. all contributions to the ith node
%from nodes j
KB = sparse(Si,Sj,Sk,NPB,NPB);
MB = sparse(Si,Sj,Mk,NPB,NPB);

ES = convhull(PB);
NES = size(ES,1);
NPS = length(unique(ES));
KS = spalloc(NPB,NPB,NPS*7);    %discrete Laplace-Beltrami
MS = spalloc(NPB,NPB,NPS*7);    %mass matrix 
% Mass matrix on reference triangle in 2D
MrefS = (ones(3)+eye(3))/24;

for i=1:NES
   e = ES(i,:);
   % vertices of the considered triangle
   P1 = PB(e(1),:);
   P2 = PB(e(2),:);
   P3 = PB(e(3),:);
   % edges of the considered triangle
   v1 = P3-P2;
   v2 = P1-P3;
   v3 = P2-P1;
   % heights of the considered triangle
   h = [v2-(v1*v2')*v1/sum(v1.^2);
        v3-(v2*v3')*v2/sum(v2.^2);
        v1-(v3*v1')*v3/sum(v3.^2)];
   % gradients of the local basis functions
   nabla = [h(1,:)/sum(h(1,:).^2);
            h(2,:)/sum(h(2,:).^2);
            h(3,:)/sum(h(3,:).^2)];
   Aloc = nabla*nabla';     %local discrete laplace-beltrami
   S = norm(cross(v1,v2));  %double area of the given triangle
   KS(e,e) = KS(e,e) + Aloc*S; %#ok
   MS(e,e) = MS(e,e) + MrefS*S; %#ok
end

KS = KS/2;

R = spalloc(NPB, NPS, NPS);
R(unique(ES), 1:NPS) = speye(NPS); % reduction matrix

MS = R'*MS*R;
KS = R'*KS*R;


end