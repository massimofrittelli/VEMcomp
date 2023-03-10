%
% Generates polyhedral mesh and matrices on the unit sphere
% How to test it: volume of the sphere = 4/3*pi*R^3
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INPUTS:
%
% - Nx: number of nodes along each direction of the bounding box
%
% OUTPUTS:
%
% - P: array of nodes
% - h: meshsize
% - K,M: stiffness and mass matrices in the bulk
% - KS,MS,CMS: stiffness, mass, and consistency matrices on the surface
% - Elements: polyhedral elements in element3d_dummy format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%function [P, h, K, M, KS, MS, CMS, boundarynode, EGamma] = tetrahedral_mesh_sphere(N)
function [P, K, ES] = tetrahedral_mesh_sphere(N)

rng(0,'twister')
rvals = linspace(-1,1,N)';
elevation = asin(rvals);
azimuth = 2*pi*linspace(0,1,N)';
radii = 3*(linspace(0,1,N)'.^(1/3));
[X,Y,Z] = sph2cart(azimuth,elevation,radii);
P = [X,Y,Z];

ES = convhull(P);
S = unique(ES);
P(S,:) = P(S,:)./vecnorm(P(S,:),2,2);
EB = delaunayTriangulation(P);
nB = size(P,1);
nS = size(S,1);
K = spalloc(nB,nB,27*nB);
for i=1:size(EB,1)
    e = EB(i,:);
    Ke = local_matrices_tet(P(e,:));
    K(e,e) = K(e,e) + Ke;
end

end