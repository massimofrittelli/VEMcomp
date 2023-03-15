function [P, h, CubicElements, NonCubicElements] = generate_mesh_3D_domain(fun, xmax, Nx, tol)
%GENERATE_MESH_3D_DOMAIN Summary of this function goes here
%   Detailed explanation goes here

Ncube = Nx^3;
hx = 2*xmax/(Nx-1); % Discretisation step along each dimension
h = hx*sqrt(3); % Meshsize

% GENERATE ONE CUBIC ELEMENT
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
x = linspace(-xmax,xmax,Nx); % Gridpoints in [-1,1]
P = zeros(Ncube,3); % Gridpoints of bounding box
for i=0:Nx-1
    for j=0:Nx-1
        for k=0:Nx-1
            P(Nx^2*i+Nx*j+k+1,:) = [x(i+1) x(j+1) x(k+1)];
        end
    end
end


% GENERATE ELEMENTS
CubicElements = [];
NonCubicElements = [];
accepted_node = false(Ncube,1);
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
            NewCubicElement = shiftElement(ESD, P(indexes(1),:));
            if is_outside(NewCubicElement, fun)
                continue
            end
            if is_inside(NewCubicElement, fun)
                % Store cubic elements that are inside domain
                accepted_node(indexes) = true(8,1);
                NewCubicElement.Pind = indexes;
                CubicElements = [CubicElements; NewCubicElement]; %#ok 
                continue
            end
            % Store non-cubic elements obtained by cutting cubic elements
            % with boundary. Such non-cunic elements are not endowed with 
            % node indexes yet.
            NewElement = cutElement(NewCubicElement, fun, tol);
            NonCubicElements = [NonCubicElements; NewElement]; %#ok 
        end
    end
end

% AFTER ELIMINATING CUBIC ELEMENTS THAT ARE OUTSIDE DOMAIN, RE-DETERMINE
% INDEXES OF NODES OF CUBIC ELEMENTS
P = P(accepted_node,:);
acceptedindexes = zeros(Ncube,1);
acceptedindexes(accepted_node,1) = linspace(1,length(P),length(P))';
for i=1:length(CubicElements)
   CubicElements(i).Pind = acceptedindexes(CubicElements(i).Pind); %#ok
end

% DETERMINE SET OF NON-REPEATED NODES UP TO SMALL TOLERANCE
for i=1:length(NonCubicElements)
   P = [P; NonCubicElements(i).P]; %#ok 
end
P = uniquetol(P,tol,'ByRows',true);

% FIX ELEMENTS BY ASSIGNING NODE INDEXES AND ELIMINATING DUPLICATE NODES UP
% TO SMALL TOLERANCE
for i=1:length(CubicElements)
   [~, ind] = ismembertol(CubicElements(i).P,P,tol,'ByRows',true);
   CubicElements(i).Pind = ind; %#ok 
   CubicElements(i).P = P(ind,:); %#ok
   for j=1:length(CubicElements(i).Faces)
       [~, ind] = ismembertol(CubicElements(i).Faces(j).P,P,tol,'ByRows',true);
       CubicElements(i).Faces(j).Pind = ind;
       CubicElements(i).Faces(j).P = P(ind,:);
   end
end

for i=1:length(NonCubicElements)
   [~, ind] = ismembertol(NonCubicElements(i).P,P,tol,'ByRows',true);
   NonCubicElements(i).Pind = ind; %#ok 
   NonCubicElements(i).P = P(ind,:); %#ok
   for j=1:length(NonCubicElements(i).Faces)
       [~, ind] = ismembertol(NonCubicElements(i).Faces(j).P,P,tol,'ByRows',true);
       NonCubicElements(i).Faces(j).Pind = ind;
       NonCubicElements(i).Faces(j).P = P(ind,:);
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
function CutElement = cutElement(CubicElement, fun, tol)
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
    if length(CutFacesP) < 4
        CutElement = [];
        return
    end
    chull = convhull(CutFacesP,'simplify',true); 
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
            CutFaces = [CutFaces; element2d_dummy(CutFacesP(chull(i,:),:), false)]; %#ok
        end
    end
    CutElement = element3d_dummy(CutFacesP, CutFaces, false);
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
    CutFace = element2d_dummy(Pnew(sort(ind),:), false);
end

% COMPUTES INTERSECTION POINT BETWEEN AN EDGE AND THE BOUNDARY OF THE
% DOMAIN
function newP = intersectEdge(P1, P2, fun)
    fun_restricted = @(alpha) fun(alpha*P1 + (1-alpha)*P2);
    alpha = fzero(fun_restricted, [0 1]);
    newP = alpha*P1 + (1-alpha)*P2;
end