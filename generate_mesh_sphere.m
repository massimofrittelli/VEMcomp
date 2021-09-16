%
% Generates polyhedral mesh on the unit sphere
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INPUTS:
%
% - Nx: number of nodes along each direction of the bounding box
%
% OUTPUTS:
%
% - P: array of nodes
% - EL: array of cuboidal elements
% - h: meshsize
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [P, EL, h] = generate_mesh_sphere(Nx)

% CREATING GRIDPOINTS OF BOUNDING BOX

hx = 1/(Nx-1); % Discretisation step along each dimension
h = hx*sqrt(3); % Meshsize
x = linspace(-1,1,Nx); % Gridpoints in [-1,1]
P = zeros(Nx^3,3); % Gridpoints of bounding box
radii = zeros(Nx^3,1); % Distances of nodes from origin (sphere centre)
for i=0:Nx-1
    for j=0:Nx-1
        for k=0:Nx-1
            P(Nx^2*i+Nx*j+k+1,:) = [x(i+1) x(j+1) x(k+1)];
            radii(Nx^2*i+Nx*j+k+1) = norm([x(i+1) x(j+1) x(k+1)]);
        end
    end
end

% CREATING CUBIC ELEMENTS AS AN INDEX ARRAY
% ELIMINATING ELEMENTS THAT ARE OUTSIDE
% RESHAPING ELEMENTS THAT INTERSECT BOUNDARY
% UPDATING DISTANCE OF NODES FROM ORIGIN

EL = zeros((Nx-1)^3,8);
outside = false((Nx-1)^3,1);
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
            if min(radii(indexes)) >= 1
                outside((Nx-1)^2*i + (Nx-1)*j + k + 1) = true;
            else
                EL((Nx-1)^2*i + (Nx-1)*j + k + 1,:) = indexes;
                if max(radii(indexes)) > 1
                    ii = indexes(radii(indexes)>1);
                    P(ii,:) = P(ii,:) ./ repmat(radii(ii),[1,3]);
                    radii(ii,:) = ones(length(ii),1);
                end
            end
        end          
    end
end
EL(outside,:) = [];

% REORDERING NODES: OUTSIDE NODE LAST

P = P(radii <= 1,:);
Ninside = size(P,1);
node_ordering(radii <= 1) = linspace(1,Ninside, Ninside)';
node_ordering(radii > 1) = linspace(1+Ninside, Nx^3, Nx^3- Ninside)';
EL = node_ordering(EL);

end