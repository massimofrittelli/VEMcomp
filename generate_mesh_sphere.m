%
% DESCRIPTION - Generates polyhedral mesh on the unit sphere throuhh the
% extrusion technique - needs fixing
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

function [P, h, SurfaceElements, BulkElements] = generate_mesh_sphere(Nx, tol)

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
for i=0:Nx-1
    for j=0:Nx-1
        for k=0:Nx-1
            P(Nx^2*i+Nx*j+k+1,:) = [x(i+1) x(j+1) x(k+1)];
        end
    end
end

% MATRIX ASSEMBLY
acceptednode = false(2*size(P,1),1);
Pnew = [];
% Twice the amount of nodes of the bounding box to allow for extrusion
BulkElements = [];
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
                BulkElements = [BulkElements; NewCubicElement]; %#ok
                
                % Extrude element, if it has an external face
                if norm(TPO) < 1 - h
                    continue
                end
                NewElements = extrude(NewCubicElement,Ncube);
                BulkElements = [BulkElements; NewElements]; %#ok
                for l=1:length(NewElements)
                    Pnew = [Pnew; NewElements(l).P]; %#ok
                end
            end
        end
    end
end

% DISCARD UNUSED NODES OF THE BOUNDING BOX
P = P(acceptednode,:);

% DETERMINE SET OF NON-REPEATED NODES UP TO SMALL TOLERANCE
P = uniquetol([P; Pnew],tol,'ByRows',true);

% FIX ELEMENTS BY ASSIGNING NODE INDEXES AND ELIMINATING DUPLICATE NODES UP
% TO SMALL TOLERANCE
SurfaceElements = [];
for i=1:length(BulkElements)
   [~, ind] = ismembertol(BulkElements(i).P,P,tol,'ByRows',true);
   BulkElements(i).Pind = ind; %#ok 
   BulkElements(i).P = P(ind,:); %#ok
   for j=1:length(BulkElements(i).Faces)
       [~, ind] = ismembertol(BulkElements(i).Faces(j).P,P,tol,'ByRows',true);
       BulkElements(i).Faces(j).Pind = ind;
       BulkElements(i).Faces(j).P = P(ind,:);
       if BulkElements(i).Faces(j).is_boundary && length(BulkElements(i).Faces(j).Pind) == 3
           SurfaceElements = [SurfaceElements; BulkElements(i).Faces(j).Pind']; %#ok
       end
   end
end

end

% Extracts largest magnitude entry from array
function y = maxabs(x)
    [m,mind] = max(abs(x));
    y = m * sign(x(mind));
end