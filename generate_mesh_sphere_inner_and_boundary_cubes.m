%
% DESCRIPTION - Generates cubic inner mesh for the unit sphere
% Notice: this function serves only as an auxiliary function for the script
% NUMPDE_MESH_GEN_STEP_2.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INPUTS:
%
% - Nx: number of nodes along each direction of the bounding box
%
% OUTPUTS:
%
% - P: array of nodes
% - h: meshsize
% - Elements: polyhedral elements in element3d_dummy format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [P, h, boundarynode, Elements] = generate_mesh_sphere_inner_and_boundary_cubes(Nx)

hx = 2/(Nx-1); % Discretisation step along each dimension
h = hx*sqrt(3); % Meshsize
Ncube = Nx^3; % Amount of nodes of the bounding box

P1S = [0 0 0; 0 1 0; 1 1 0; 1 0 0]*hx; % bottom face
P2S = [0 0 1; 0 1 1; 1 1 1; 1 0 1]*hx; % top face
P3S = [0 0 0; 0 1 0; 0 1 1; 0 0 1]*hx; % back face
P4S = [1 0 0; 1 1 0; 1 1 1; 1 0 1]*hx; % front face
P5S = [0 0 0; 1 0 0; 1 0 1; 0 0 1]*hx; % left face
P6S = [0 1 0; 1 1 0; 1 1 1; 0 1 1]*hx; % right face
PS = unique([P1S; P2S; P3S; P4S; P5S; P6S],'rows');

E1S = element2d(P1S, true);
E2S = element2d(P2S, true);
E3S = element2d(P3S, true);
E4S = element2d(P4S, true);
E5S = element2d(P5S, true);
E6S = element2d(P6S, true);
ESD = element3d(PS, [E1S;E2S;E3S;E4S;E5S;E6S], true);

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

% GENERATE ELEMENTS
acceptednode = false(2*size(P,1),1);
boundarynode = false(2*size(P,1),1);
newP = [P; zeros(size(P))];
% Twice the amount of nodes of the bounding box to allow for extrusion
Elements = [];
for i=0:Nx-2 % For each element of the bounding box
    for j=0:Nx-2
        for k=0:Nx-2
            % Inner test point to determine if cubic element is fully
            % outside the shpere
             TPI = [minabs([x(i+1) x(i+2)]) minabs([x(j+1) x(j+2)]) minabs([x(k+1) x(k+2)])];

            if norm(TPI) < 1
                % Element not fully outside, we keep it
                % Otherwise we discard it
                indexes = [Nx^2*i+Nx*j+k+1
                       Nx^2*i+Nx*j+k+2
                       Nx^2*i+Nx*(j+1)+k+1
                       Nx^2*i+Nx*(j+1)+k+2
                       Nx^2*(i+1)+Nx*j+k+1
                       Nx^2*(i+1)+Nx*j+k+2
                       Nx^2*(i+1)+Nx*(j+1)+k+1
                       Nx^2*(i+1)+Nx*(j+1)+k+2];
                acceptednode(indexes,1) = true(8,1);                   
                NewCubicElement = shiftElement(ESD, P(indexes(1),:));
                Elements = [Elements; NewCubicElement]; %#ok 
            end
        end
    end
end


P = newP(acceptednode,:);
boundarynode = find(boundarynode(acceptednode));


end

function y = minabs(x)
    [m,mind] = min(abs(x));
    y = m * sign(x(mind));
end