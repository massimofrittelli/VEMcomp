function [ A,M ] = matrices( P,E )

% DESCRIPTION - Generates SFEM matrices for an arbitrary surface mesh.

% INPUTS:
%  - P is a NPx3 array containing the coordinates of the nodes
%  - E is a NEx3 array in which each row contains the indexes of the nodes
%    of a triangular element (common representation of elements)

% OUTPUTS:
%  - A is the NPxNP stiffness matrix for the Laplace-Beltrami operator
%  - M is the mass matrix

NE = size(E,1);
NP = size(P,1);
A = spalloc(NP,NP,NP*7);    %stiffness matrix for Laplace-Beltrami
M = spalloc(NP,NP,NP*7);    %mass matrix       
refM = (ones(3)+eye(3))/24; %mass matrix on reference element

for i=1:NE
   e = E(i,:);
   % vertices of the considered triangle
   P1 = P(e(1),:);
   P2 = P(e(2),:);
   P3 = P(e(3),:);
   % edges of the considered triangle
   v1 = P3-P2;
   v2 = P1-P3;
   v3 = P2-P1;
   % heights of the considered triangle
   h = [v2-(v1*v2')*v1/norm(v1)^2;
        v3-(v2*v3')*v2/norm(v2)^2;
        v1-(v3*v1')*v3/norm(v3)^2];
   % gradients of the local basis functions
   nabla = [h(1,:)/norm(h(1,:))^2;
            h(2,:)/norm(h(2,:))^2;
            h(3,:)/norm(h(3,:))^2];
   Aloc = nabla*nabla';     %local discrete laplace-beltrami
   S = norm(cross(v1,v2));  %double area of the given triangle
   A(e,e) = A(e,e) + Aloc*S;
   M(e,e) = M(e,e) + refM*S;
end

A = A/2;

end