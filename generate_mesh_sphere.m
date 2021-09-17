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

function [P, h, K, M, radii, Elements] = generate_mesh_sphere(Nx)

hx = 2/(Nx-1); % Discretisation step along each dimension
h = hx*sqrt(3); % Meshsize

% COMPUTE MATRICES ON REFERENCE CUBE
K = spalloc(Nx^3,Nx^3,57*Nx^3);
M = spalloc(Nx^3,Nx^3,57*Nx^3);

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

% MATRIX ASSEMBLY
newP = P;
Elements = [];
for i=0:Nx-2 % For each element of the bounding box
    for j=0:Nx-2
        for k=0:Nx-2
            % Create test points to determine if cubic element is fully
            % outside, fully inside or intersects the boundary.
            % Outer test point
            TPO = [max(abs([x(i+1) x(i+2)])) max(abs([x(j+1) x(j+2)])) max(abs([x(k+1) x(k+2)]))];
            % Inner test point
            TPI = [min(abs([x(i+1) x(i+2)])) min(abs([x(j+1) x(j+2)])) min(abs([x(k+1) x(k+2)]))];
            if norm(TPI) >= 1
                % Element fully outside, discard
                continue
            else
                indexes = [Nx^2*i+Nx*j+k+1
                       Nx^2*i+Nx*j+k+2
                       Nx^2*i+Nx*(j+1)+k+1
                       Nx^2*i+Nx*(j+1)+k+2
                       Nx^2*(i+1)+Nx*j+k+1
                       Nx^2*(i+1)+Nx*j+k+2
                       Nx^2*(i+1)+Nx*(j+1)+k+1
                       Nx^2*(i+1)+Nx*(j+1)+k+2];
                if norm(TPO) <= 1
                    % Element fully inside, is a cube
                    M(indexes, indexes) = M(indexes, indexes) + MS; %#ok
                    K(indexes, indexes) = K(indexes, indexes) + KS; %#ok
                    
                    Elements = [Elements; shiftElement(ESD, P(indexes(1),:))]; %#ok
                else
                    % Element intersects the boundary, cut cube                    
                    
                    P1 = P(indexes([1 3 7 5]),:);
                    P2 = P(indexes([2 4 8 6]),:);
                    P3 = P(indexes([1 3 4 2]),:);
                    P4 = P(indexes([5 7 8 6]),:);
                    P5 = P(indexes([1 5 6 2]),:);
                    P6 = P(indexes([3 7 8 4]),:);
                    
                    E1D = element2ddummy(P1, true);
                    E2D = element2ddummy(P2, true);
                    E3D = element2ddummy(P3, true);
                    E4D = element2ddummy(P4, true);
                    E5D = element2ddummy(P5, true);
                    E6D = element2ddummy(P6, true);

                    ED = element3ddummy(P(indexes,:), [E1D;E2D;E3D;E4D;E5D;E6D], true);
                    EC = cut(ED);
                   
                    E = dummy2element(EC);
                    Elements = [Elements; EC]; %#ok
                                        
                    M(indexes, indexes) = M(indexes, indexes) + E.M; %#ok
                    K(indexes, indexes) = K(indexes, indexes) + E.K; %#ok
                    
                    newP(indexes,:) = EC.P;
                end
            end
        end
    end
end

P = newP;



end