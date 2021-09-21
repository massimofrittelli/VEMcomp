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

function [P, h, K, M, Elements] = generate_mesh_sphere(Nx)

hx = 2/(Nx-1); % Discretisation step along each dimension
h = hx*sqrt(3); % Meshsize
Ncube = Nx^3; % Amount of nodes of the bounding box

% COMPUTE MATRICES ON REFERENCE CUBE
K = spalloc(2*Ncube,2*Ncube,57*Ncube);
M = spalloc(2*Ncube,2*Ncube,57*Ncube);
% Twice the amount of nodes of the bounding box to allow for extrusion

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
ES = dummy2element(ESD);

KS = ES.K;
MS = ES.M;

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
acceptednode = false(2*size(P,1),1);
newP = [P; zeros(size(P))];
% Twice the amount of nodes of the bounding box to allow for extrusion
Elements = [];
for i=0:Nx-2 % For each element of the bounding box
    for j=0:Nx-2
        for k=0:Nx-2
            % Outer test point to determine if cubic element is fully
            % inside the shpere
            TPO = [maxabs([x(i+1) x(i+2)]) maxabs([x(j+1) x(j+2)]) maxabs([x(k+1) x(k+2)])];
            if norm(TPO) <= 1
                % Element fully inside, we keep it
                % Otherwise we discard it
                % TODO: use dummy 3d element to save operations
                indexes = [Nx^2*i+Nx*j+k+1
                       Nx^2*i+Nx*j+k+2
                       Nx^2*i+Nx*(j+1)+k+1
                       Nx^2*i+Nx*(j+1)+k+2
                       Nx^2*(i+1)+Nx*j+k+1
                       Nx^2*(i+1)+Nx*j+k+2
                       Nx^2*(i+1)+Nx*(j+1)+k+1
                       Nx^2*(i+1)+Nx*(j+1)+k+2];
                acceptednode(indexes,1) = true(8,1);
                M(indexes, indexes) = M(indexes, indexes) + MS; %#ok
                K(indexes, indexes) = K(indexes, indexes) + KS; %#ok
                    
                NewCubicElement = shiftElement(ESD, P(indexes(1),:));
                Elements = [Elements; NewCubicElement]; %#ok
                
                % Extrude element, if it has an external face
                extrusion_directions = find(vecnorm(repmat(TPO',1,3) + eye(3)*hx) > 1 | vecnorm(repmat(TPO',1,3) - eye(3)*hx) > 1);
                if isempty(extrusion_directions)
                   continue 
                end
                extrusion_verses = sign(TPO(extrusion_directions));
                NewElements = extrude(NewCubicElement,indexes,extrusion_directions,extrusion_verses,Ncube);
                Elements = [Elements; NewElements]; %#ok
                for l=1:length(NewElements)
                    eind = NewElements(l).Pind;
                    E = dummy2element(NewElements(l));
                    acceptednode(eind, 1) = true(8,1);
                    newP(eind(5:8),:) = NewElements(l).P(5:8,:);
                    M(eind, eind) = M(eind, eind) + E.M; %#ok
                    K(eind, eind) = K(eind, eind) + E.K; %#ok
%                     if abs(sum(sum(E.M)) - newP(eind,:)'*E.K*newP(eind,:)) > 1e-14
%                         error('Error!')
%                     end
                end
            end
        end
    end
end

P = newP(acceptednode,:);
M = M(acceptednode, acceptednode);
K = K(acceptednode, acceptednode);


end

% Extracts largest magnitude element from array
function y = maxabs(x)
    [m,mind] = max(abs(x));
    y = m * sign(x(mind));
end