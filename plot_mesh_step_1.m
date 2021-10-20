%
% Generates polyhedral mesh on the unit sphere
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
% - K,M: stiffness and mass matrices
% - Elements: polyhedral elements in element3ddummy format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [P, h, K, M, KS, MS, boundarynode, EGamma, Elements] = plot_mesh_step_1(Nx)

hx = 2/(Nx-1); % Discretisation step along each dimension
h = hx*sqrt(3); % Meshsize
Nsurf = 6*Nx^2-12*Nx-8; % Amount of nodes on the boundary of the bounding box
Ncube = Nx^3; % Amount of nodes of the bounding box

% COMPUTE MATRICES ON REFERENCE CUBE
K = spalloc(2*Ncube,2*Ncube,57*Ncube); % Stiffness matrix in the bulk
M = spalloc(2*Ncube,2*Ncube,57*Ncube); % Mass matrix in the bulk
% Twice the amount of nodes of the bounding box to allow for extrusion

KS = spalloc(2*Ncube,2*Ncube,9*Nsurf); % Stiffness matrix on the surf
MS = spalloc(2*Ncube,2*Ncube,9*Nsurf); % Mass matrix on the surf

P1S = [0 0 0; 0 1 0; 1 1 0; 1 0 0]*hx; % bottom face
P2S = [0 0 1; 0 1 1; 1 1 1; 1 0 1]*hx; % top face
P3S = [0 0 0; 0 1 0; 0 1 1; 0 0 1]*hx; % back face
P4S = [1 0 0; 1 1 0; 1 1 1; 1 0 1]*hx; % front face
P5S = [0 0 0; 1 0 0; 1 0 1; 0 0 1]*hx; % left face
P6S = [0 1 0; 1 1 0; 1 1 1; 0 1 1]*hx; % right face

E1S = element2ddummy(P1S, true);
E2S = element2ddummy(P2S, true);
E3S = element2ddummy(P3S, true);
E4S = element2ddummy(P4S, true);
E5S = element2ddummy(P5S, true);
E6S = element2ddummy(P6S, true);
PS = unique([P1S; P2S; P3S; P4S; P5S; P6S],'rows');
ESD = element3ddummy(PS, [E1S;E2S;E3S;E4S;E5S;E6S], true);
EC = dummy2element(ESD);

KC = EC.K;
MC = EC.M;

% CREATING GRIDPOINTS OF BOUNDING BOX

x = linspace(-1,1,Nx); % Gridpoints in [-1,1]
P = zeros(Ncube,3); % Gridpoints of bounding box
radii = zeros(Ncube,1); % Distances of nodes from origin (sphere centre)
for i=0:Nx-1
    for j=0:Nx-1
        for k=0:Nx-1
            P(Nx^2*i+Nx*j+k+1,:) = [x(i+1) x(j+1) x(k+1)];
            radii(Nx^2*i+Nx*j+k+1) = norm([x(i+1) x(j+1) x(k+1)]);
        end
    end
end

% MATRIX ASSEMBLY
boundarynode = false(2*size(P,1),1);
EGamma = [];
% Twice the amount of nodes of the bounding box to allow for extrusion
Elements = [];
for i=0:Nx-2 % For each element of the bounding box
    for j=0:Nx-2
        for k=0:Nx-2
                indexes = [Nx^2*i+Nx*j+k+1
                       Nx^2*i+Nx*j+k+2
                       Nx^2*i+Nx*(j+1)+k+1
                       Nx^2*i+Nx*(j+1)+k+2
                       Nx^2*(i+1)+Nx*j+k+1
                       Nx^2*(i+1)+Nx*j+k+2
                       Nx^2*(i+1)+Nx*(j+1)+k+1
                       Nx^2*(i+1)+Nx*(j+1)+k+2];
                M(indexes, indexes) = M(indexes, indexes) + MC; %#ok
                K(indexes, indexes) = K(indexes, indexes) + KC; %#ok
                    
                NewCubicElement = shiftElement(ESD, P(indexes(1),:));
                Elements = [Elements; NewCubicElement]; %#ok 
        end
    end
end


MS = MS(boundarynode, boundarynode);
KS = KS(boundarynode, boundarynode);


end