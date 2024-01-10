function [P, h, BulkElements, SurfElements, ElementsPlot] = generate_mesh_3d(fun, Q, Nx, tol, xcut)
%GENERATE_MESH_3D_DOMAIN Summary of this function goes here
%   Detailed explanation goes here

hx_requested = (Q(1,2)-Q(1,1))/(Nx-1); % Requested discretisation step along x
hy_requested = (Q(2,2)-Q(2,1))/(Nx-1); % Requested discretisation step along y
hz_requested = (Q(3,2)-Q(3,1))/(Nx-1); % Requested discretisation step along z

h_mono = min([hx_requested,hy_requested,hz_requested]); % Actual discretization step along each dimension
h = h_mono*sqrt(3); % Meshsize

Nx = ceil((Q(1,2)-Q(1,1))/h_mono)+1; % Corrected number of discretization nodes along x
Ny = ceil((Q(2,2)-Q(2,1))/h_mono)+1; % Corrected number of discretization nodes along y
Nz = ceil((Q(3,2)-Q(3,1))/h_mono)+1; % Corrected number of discretization nodes along z
N = Nx*Ny*Nz; % Number of nodes of bounding box

Q(1,:) = Q(1,:) + (h_mono*(Nx-1) - (Q(1,2)-Q(1,1)))/2*[-1,1]; % Corrected range of bounding box along x
Q(2,:) = Q(2,:) + (h_mono*(Ny-1) - (Q(2,2)-Q(2,1)))/2*[-1,1]; % Corrected range of bounding box along y
Q(3,:) = Q(3,:) + (h_mono*(Nz-1) - (Q(3,2)-Q(3,1)))/2*[-1,1]; % Corrected range of bounding box along z

i_cut = min(max(round((xcut-Q(1,1))/h_mono),0),Nx);

% GENERATE ONE CUBIC ELEMENT
P1S = [0 0 0; 0 1 0; 1 1 0; 1 0 0]*h_mono; % bottom face
P2S = [0 0 1; 1 0 1; 1 1 1; 0 1 1]*h_mono; % top face
P3S = [0 0 0; 0 0 1; 0 1 1; 0 1 0]*h_mono; % back face
P4S = [1 0 0; 1 1 0; 1 1 1; 1 0 1]*h_mono; % front face
P5S = [0 0 0; 1 0 0; 1 0 1; 0 0 1]*h_mono; % left face
P6S = [0 1 0; 0 1 1; 1 1 1; 1 1 0]*h_mono; % right face

E1S = element2d(P1S, true);
E2S = element2d(P2S, true);
E3S = element2d(P3S, true);
E4S = element2d(P4S, true);
E5S = element2d(P5S, true);
E6S = element2d(P6S, true);
PS = unique([P1S; P2S; P3S; P4S; P5S; P6S],'rows');

[~, p1, q1] = intersect(PS, P1S, 'rows', 'stable');
[~, p2, q2] = intersect(PS, P2S, 'rows', 'stable');
[~, p3, q3] = intersect(PS, P3S, 'rows', 'stable');
[~, p4, q4] = intersect(PS, P4S, 'rows', 'stable');
[~, p5, q5] = intersect(PS, P5S, 'rows', 'stable');
[~, p6, q6] = intersect(PS, P6S, 'rows', 'stable');

setPind(E1S, p1(q1));
setPind(E2S, p2(q2));
setPind(E3S, p3(q3));
setPind(E4S, p4(q4));
setPind(E5S, p5(q5));
setPind(E6S, p6(q6));

ESD = element3d(PS, [E1S;E2S;E3S;E4S;E5S;E6S], true, (1:8)');

% CREATING GRIDPOINTS OF BOUNDING BOX
x = linspace(Q(1,1),Q(1,2),Nx); % Gridpoints in [-1,1]
y = linspace(Q(2,1),Q(2,2),Ny); % Gridpoints in [-1,1]
z = linspace(Q(3,1),Q(3,2),Nz); % Gridpoints in [-1,1]
P = zeros(N,3); % Gridpoints of bounding box
for i=0:Nx-1
    for j=0:Ny-1
        for k=0:Nz-1
            P(Ny*Nz*i+Nz*j+k+1,:) = [x(i+1) y(j+1) z(k+1)];
        end
    end
end


% GENERATE ELEMENTS
BulkElements = [];
accepted_node = false(N,1);
newP = [];
for i=0:Nx-2 % For each element of the bounding box
    for j=0:Ny-2
        for k=0:Nz-2
            NewCubicElement = shiftElement(ESD, [x(i+1) y(j+1) z(k+1)]);
            if is_outside(NewCubicElement, fun)
                continue
            end
            if is_inside(NewCubicElement, fun)
                indexes = [Nz*Ny*i+Nz*j+k+1
                           Nz*Ny*i+Nz*j+k+2
                           Nz*Ny*i+Nz*(j+1)+k+1
                           Nz*Ny*i+Nz*(j+1)+k+2
                           Nz*Ny*(i+1)+Nz*j+k+1
                           Nz*Ny*(i+1)+Nz*j+k+2
                           Nz*Ny*(i+1)+Nz*(j+1)+k+1
                           Nz*Ny*(i+1)+Nz*(j+1)+k+2];
                % Store cubic elements that are inside domain
                accepted_node(indexes) = true(8,1);
                if i == i_cut
                   for fc =1:6
                      if max(NewCubicElement.Faces(fc).P(:,1)) == x(1+i_cut)
                         NewCubicElement.Faces(fc).to_plot = true;
                      end
                   end
                end
                BulkElements = [BulkElements; NewCubicElement]; %#ok
                continue
            end
            % Store non-cubic elements obtained by cutting cubic elements
            % with boundary. Such non-cubic elements are not endowed with 
            % node indexes yet.
            NewElement = cutElement(NewCubicElement, fun, tol);
            if not(isempty(NewElement))
                if i >= i_cut
                    for fc = 1:NewElement.NFaces
                       if NewElement.Faces(fc).is_boundary
                           NewElement.Faces(fc).to_plot = true;
                       end
                    end
                end
                if i == i_cut
                    for fc = 1:NewElement.NFaces
                       if (max(NewElement.Faces(fc).P(:,1)) - x(1+i_cut)) < tol
                           NewElement.Faces(fc).to_plot = true;
                       end
                    end
                end
                BulkElements = [BulkElements; NewElement]; %#ok
                newP = [newP; NewElement.P];  %#ok
            end
        end
    end
end

% DETERMINE SET OF NON-REPEATED NODES UP TO SMALL TOLERANCE
P = uniquetol([P(accepted_node,:); newP],tol,'ByRows',true);

ElementsPlot = [];
% FIX ELEMENTS BY ASSIGNING NODE INDEXES AND ELIMINATING DUPLICATE NODES UP
% TO SMALL TOLERANCE
SurfElements = [];
for i=1:length(BulkElements)
   [~, ind] = ismembertol(BulkElements(i).P,P,tol,'ByRows',true);
   setPind(BulkElements(i), ind);
   setP(BulkElements(i), P(ind,:));
   for j=1:length(BulkElements(i).Faces)
       [~, ind] = ismembertol(BulkElements(i).Faces(j).P,P,tol,'ByRows',true);
       setPind(BulkElements(i).Faces(j), ind);
       setP(BulkElements(i).Faces(j), P(ind,:));
       if BulkElements(i).Faces(j).is_boundary
           SurfElements = [SurfElements; ind']; %#ok
       end
       if BulkElements(i).Faces(j).to_plot
           ElementsPlot = [ElementsPlot; BulkElements(i).Faces(j)]; %#ok
       end
   end
end



end

% DETERMINE IF ELEMENT IS OUTSIDE DOMAIN (POSSIBLY TOUCHING BOUNDARY)
function outside = is_outside(Element, fun)
    outside = true;
    for i=1:length(Element.P)
        if fun(Element.P(i,:)) < 0
            outside = false;
            break
        end
    end
end

% DETERMINE IF ELEMENT IS INSIDE DOMAIN (POSSIBLY TOUCHING BOUNDARY)
function inside = is_inside(Element, fun)
    inside = true;
    for i=1:length(Element.P)
        if fun(Element.P(i,:)) > 0
            inside = false;
            break
        end
    end
end

% CUTS A CUBIC ELEMENT BY BOUNDARY OF DOMAIN
function [CutElement] = cutElement(CubicElement, fun, tol)
    CutFaces = [];
    CutFacesP = [];
    for i=1:6
        if is_inside(CubicElement.Faces(i), fun)
           CutFaces = [CutFaces; CubicElement.Faces(i)]; %#ok
           CutFacesP = [CutFacesP; CubicElement.Faces(i).P]; %#ok
           continue
        end
        if not(is_outside(CubicElement.Faces(i), fun))
           newCutFace = cutFace(CubicElement.Faces(i), fun, tol);
           if not(isempty(newCutFace))
                CutFaces = [CutFaces; newCutFace]; %#ok
                CutFacesP = [CutFacesP; newCutFace.P]; %#ok
           end
        end
    end
    CutFacesP = uniquetol(CutFacesP,tol,'ByRows',true);
    if length(CutFacesP) < 4 % Element is flat -> discard
        CutElement = [];
        return
    end
    try
        chull = convhull(CutFacesP,'simplify',true); 
    catch % Element has zero volume -> discard
        warning('Almost-zero volume element discarded -> tiny dent in the mesh')
        CutElement = [];
        return
    end
    for i=1:size(chull,1)
        contained = false;
        for j=1:length(CutFaces)
            ind = ismembertol(CutFacesP(chull(i,:),:), CutFaces(j).P,tol,'ByRows',true);
            if nnz(ind) == 3
                contained = true;
                break
            end
        end
        if not(contained)
            newFaceP = CutFacesP(chull(i,:),:);
            CutFaces = [CutFaces; element2d(newFaceP, false, true, [], mean(newFaceP,1))]; %#ok
        end
    end
    CutElement = element3d(CutFacesP, CutFaces, false, [], mean(CutFacesP,1));
end

% CUTS A FACE BY BOUNDARY OF DOMAIN
function CutFace = cutFace(Face, fun, tol)
    P = Face.P;
    inside_or_boundary = fun(P) <= 0;
    Pnew = [];
    for i=1:length(P)
       j = rem(i, length(P)) + 1;
       if inside_or_boundary(i)
           Pnew = [Pnew; P(i,:)]; %#ok
       end
       if sum(inside_or_boundary([i j])) == 1
           Pnew = [Pnew; intersectEdge(P(i,:), P(j,:), fun)]; %#ok
       end
    end
    [~, ind] = uniquetol(Pnew,tol,'Byrows',true);
    if length(ind) < 3
        CutFace = [];
        return
    end
    cutFaceP = Pnew(sort(ind),:);
    CutFace = element2d(cutFaceP, false, false, [], mean(cutFaceP,1));
end

% COMPUTES INTERSECTION POINT BETWEEN AN EDGE AND THE BOUNDARY OF THE
% DOMAIN
function newP = intersectEdge(P1, P2, fun)
    fun_restricted = @(alpha) fun(alpha*P1 + (1-alpha)*P2);
    alpha = fzero(fun_restricted, [0 1]);
    newP = alpha*P1 + (1-alpha)*P2;
end