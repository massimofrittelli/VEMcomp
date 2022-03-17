%
% DESCRIPTION - Generates polyhedral mesh and VEM matrices on the unit sphere
% Hint for testing: volume of the sphere = 4/3*pi*R^3
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
% - boundarynode: boolean array that identifies nodes that are on the surface
% - EGamma: all elements of the surface
% - ElementsCut: elements of the bulk. Not all of them, just the ones in
%   proximity of the cut (sphere is cut to see inside)
% - EGammaCut: elements of the surface. Not all of them, just the one that
%   survive after the sphere is cut.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [P, h, K, M, KS, MS, CMS, boundarynode, EGamma, ElementsCut, EGammaCut] = generate_mesh_sphere(Nx)

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
CMS = spalloc(2*Ncube,2*Ncube,9*Nsurf); % Consistency matrix on the surf

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
acceptednode = false(2*size(P,1),1);
boundarynode = false(2*size(P,1),1);
newP = [P; zeros(size(P))];
EGamma = [];
EGammaCut = [];
% Twice the amount of nodes of the bounding box to allow for extrusion
ElementsCut = [];
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
                M(indexes, indexes) = M(indexes, indexes) + MC; %#ok
                K(indexes, indexes) = K(indexes, indexes) + KC; %#ok
                    
                NewCubicElement = shiftElement(ESD, P(indexes(1),:));
                if i == ceil((Nx-2)/2) || j == ceil((Nx-2)/2)
                    NewCubicElement.Pind = indexes;
                    ElementsCut = [ElementsCut; NewCubicElement]; %#ok
                end
                
                % Extrude element, if it has an external face
                if norm(TPO) < 1 - h
                    continue
                end
                NewCubicElement.Pind = indexes;
                NewElements = extrude(NewCubicElement,Ncube);
                if i == ceil((Nx-2)/2) || j == ceil((Nx-2)/2)
                    ElementsCut = [ElementsCut; NewElements]; %#ok
                end
                for l=1:length(NewElements)
                    eind = NewElements(l).Pind;
                    E = dummy2element(NewElements(l));
                    acceptednode(eind, 1) = true(8,1);
                    boundarynode(eind(5:8),1) = true(4,1);
                    EGamma = [EGamma; eind(5:7)'; eind([5 7:8])']; %#ok
                    if i >= ceil((Nx-2)/2) || j >= ceil((Nx-2)/2)
                        EGammaCut = [EGammaCut; eind(5:7)'; eind([5 7:8])']; %#ok
                    end
                    newP(eind(5:8),:) = NewElements(l).P(5:8,:);
                    M(eind, eind) = M(eind, eind) + E.M; %#ok
                    K(eind, eind) = K(eind, eind) + E.K; %#ok
%                     if abs(sum(sum(E.M)) - newP(eind,:)'*E.K*newP(eind,:)) > 1e-14
%                         error('Error!')
%                     end
                    MS(eind(5:7), eind(5:7)) = MS(eind(5:7), eind(5:7)) + E.Faces(2).M; %#ok
                    MS(eind([5 7:8]), eind([5 7:8])) = MS(eind([5 7:8]), eind([5 7:8])) + E.Faces(3).M; %#ok
                    KS(eind(5:7), eind(5:7)) = KS(eind(5:7), eind(5:7)) + E.Faces(2).K; %#ok
                    KS(eind([5 7:8]), eind([5 7:8])) = KS(eind([5 7:8]), eind([5 7:8])) + E.Faces(3).K; %#ok
                    CMS(eind(5:7), eind(5:7)) = CMS(eind(5:7), eind(5:7)) + E.Faces(2).CM; %#ok
                    CMS(eind([5 7:8]), eind([5 7:8])) = CMS(eind([5 7:8]), eind([5 7:8])) + E.Faces(3).CM; %#ok
                end
            end
        end
    end
end


P = newP(acceptednode,:);
M = M(acceptednode, acceptednode);
K = K(acceptednode, acceptednode);
MS = MS(boundarynode, boundarynode);
KS = KS(boundarynode, boundarynode);
CMS = CMS(boundarynode, boundarynode);
boundarynode = find(boundarynode(acceptednode));
acceptedindexes = zeros(2*Ncube,1);
acceptedindexes(acceptednode,1) = linspace(1,length(P),length(P))';
EGamma = acceptedindexes(EGamma);
EGammaCut = acceptedindexes(EGammaCut);

for i=1:length(ElementsCut)
   ElementsCut(i).Pind = acceptedindexes(ElementsCut(i).Pind); %#ok
   for j=1:length(ElementsCut(i).Faces)
       ElementsCut(i).Faces(j).Pind = acceptedindexes(ElementsCut(i).Faces(j).Pind);
   end
end


end

% Extracts largest magnitude entry from array
function y = maxabs(x)
    [m,mind] = max(abs(x));
    y = m * sign(x(mind));
end