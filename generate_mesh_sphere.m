%
% DESCRIPTION - Generates polyhedral mesh on the unit sphere
% Notice: this function serves only as an auxiliary function for the
% scripts MESH_PLOTTER_STEP_3 and MESH_PLOTTER_STEP_4
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
% - EGamma: all elements of the surface
% - ElementsCut: elements of the bulk. Not all of them, just the ones in
%   proximity of the cut (sphere is cut to see inside)
% - EGammaCut: elements of the surface. Not all of them, just the one that
%   survive after the sphere is cut.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [P, h, boundarynode, EGamma, Elements, EGammaCut, ElementsCut] = generate_mesh_sphere(Nx)

hx = 2/(Nx-1); % Discretisation step along each dimension
h = hx*sqrt(3); % Meshsize
% Nsurf = 6*Nx^2-12*Nx-8; % Amount of nodes on the boundary of the bounding box
Ncube = Nx^3; % Amount of nodes of the bounding box

P1S = [0 0 0; 0 1 0; 1 1 0; 1 0 0]*hx; % bottom face
P2S = [0 0 1; 1 0 1; 1 1 1; 0 1 1]*hx; % top face
P3S = [0 0 0; 0 0 1; 0 1 1; 0 1 0]*hx; % back face
P4S = [1 0 0; 1 1 0; 1 1 1; 1 0 1]*hx; % front face
P5S = [0 0 0; 1 0 0; 1 0 1; 0 0 1]*hx; % left face
P6S = [0 1 0; 0 1 1; 1 1 1; 1 1 0]*hx; % right face

E1S = element2d_dummy(P1S, true);
E2S = element2d_dummy(P2S, true);
E3S = element2d_dummy(P3S, true);
E4S = element2d_dummy(P4S, true);
E5S = element2d_dummy(P5S, true);
E6S = element2d_dummy(P6S, true);
PS = unique([P1S; P2S; P3S; P4S; P5S; P6S],'rows');

[~, p1, q1] = intersect(PS, P1S, 'rows', 'stable');
[~, p2, q2] = intersect(PS, P2S, 'rows', 'stable');
[~, p3, q3] = intersect(PS, P3S, 'rows', 'stable');
[~, p4, q4] = intersect(PS, P4S, 'rows', 'stable');
[~, p5, q5] = intersect(PS, P5S, 'rows', 'stable');
[~, p6, q6] = intersect(PS, P6S, 'rows', 'stable');

E1S.Pind = p1(q1);
E2S.Pind = p2(q2);
E3S.Pind = p3(q3);
E4S.Pind = p4(q4);
E5S.Pind = p5(q5);
E6S.Pind = p6(q6);

ESD = element3d_dummy(PS, [E1S;E2S;E3S;E4S;E5S;E6S], true);
ESD.Pind = (1:8)';

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
% Twice the amount of nodes of the bounding box to allow for extrusion
Elements = [];
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
                NewCubicElement.Pind = indexes;
                Elements = [Elements; NewCubicElement]; %#ok
                if i == ceil((Nx-2)/2) || j == ceil((Nx-2)/2)
                    ElementsCut = [ElementsCut; NewCubicElement]; %#ok
                end
                
                % Extrude element, if it has an external face
                if norm(TPO) < 1 - h
                    continue
                end
                NewElements = extrude(NewCubicElement,Ncube);
                Elements = [Elements; NewElements]; %#ok
                if i == ceil((Nx-2)/2) || j == ceil((Nx-2)/2)
                    ElementsCut = [ElementsCut; NewElements]; %#ok
                end
                for l=1:length(NewElements)
                    Element = NewElements(l);
                    eind = Element.Pind;
                    eind_boundary = get_P_indexes_boundary(Element);
                    eind_boundary_1 = Element.Faces(2).Pind;
                    eind_boundary_2 = Element.Faces(3).Pind;
                    acceptednode(eind, 1) = true(length(eind),1);
                    boundarynode(eind_boundary,1) = true(4,1);
                    EGamma = [EGamma; eind_boundary_1'; eind_boundary_2']; %#ok
                    if i >= ceil((Nx-2)/2) || j >= ceil((Nx-2)/2)
                        EGammaCut = [EGammaCut; eind_boundary_1'; eind_boundary_2']; %#ok
                    end
                    [~,id1, id2] = intersect(eind,eind_boundary,'stable');
                    newP(eind_boundary,:) = Element.P(id1(id2),:);
                end
            end
        end
    end
end


P = newP(acceptednode,:);
boundarynode = find(boundarynode(acceptednode));
acceptedindexes = zeros(2*Ncube,1);
acceptedindexes(acceptednode,1) = linspace(1,length(P),length(P))';
EGamma = acceptedindexes(EGamma);
EGammaCut = acceptedindexes(EGammaCut);

for i=1:length(Elements)
   Elements(i).Pind = acceptedindexes(Elements(i).Pind); %#ok
   for j=1:length(Elements(i).Faces)
       Elements(i).Faces(j).Pind = acceptedindexes(Elements(i).Faces(j).Pind);
   end
end

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