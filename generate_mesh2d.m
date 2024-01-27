function [P, h, BulkElements, SurfElements] = generate_mesh2d(fun, Q, Nx, tol)
%GENERATE_MESH_FLAT_DOMAIN Summary of this function goes here
%   Detailed explanation goes here

hx_requested = (Q(1,2)-Q(1,1))/(Nx-1); % Requested discretisation step along x
hy_requested = (Q(2,2)-Q(2,1))/(Nx-1); % Requested discretisation step along y

h_mono = min([hx_requested,hy_requested]); % Actual discretization step along each dimension
h = h_mono*sqrt(2); % Meshsize

Nx = ceil((Q(1,2)-Q(1,1))/h_mono)+1; % Corrected number of discretization nodes along x
Ny = ceil((Q(2,2)-Q(2,1))/h_mono)+1; % Corrected number of discretization nodes along y
Nrect = Nx*Ny; % Number of nodes of bounding box

Q(1,:) = Q(1,:) + (h_mono*(Nx-1) - (Q(1,2)-Q(1,1)))/2*[-1,1]; % Corrected range of bounding box along x
Q(2,:) = Q(2,:) + (h_mono*(Ny-1) - (Q(2,2)-Q(2,1)))/2*[-1,1]; % Corrected range of bounding box along y

% GENERATE ONE SQUARE ELEMENT
PS = [0 0 0; 0 1 0; 1 1 0; 1 0 0]*h_mono;
ES = element2d(PS, true);

% CREATING GRIDPOINTS OF BOUNDING BOX
x = linspace(Q(1,1),Q(1,2),Nx); % Gridpoints along x
y = linspace(Q(2,1),Q(2,2),Ny); % Gridpoints along y
P = zeros(Nrect,3); % Gridpoints of bounding box
for i=0:Nx-1
    for j=0:Nx-1  
    	P(Nx*i+j+1,:) = [x(i+1) y(j+1) 0];
    end
end

% GENERATE ELEMENTS
SquareElements = [];
NonSquareElements = [];
SE = cell(0,1);
accepted_node = false(Nrect,1);
for i=0:Nx-2 % For each element of the bounding box
    for j=0:Ny-2
        indexes = [Ny*i+j+1
                   Ny*i+j+2
                   Ny*(i+1)+j+2
                   Ny*(i+1)+j+1];                  
        NewSquareElement = shiftElement(ES, P(indexes(1),:));
        if is_outside(NewSquareElement, fun)
            continue
        end
        if is_inside(NewSquareElement, fun)
            % Store square elements that are inside domain
            accepted_node(indexes) = true(4,1);
            setPind(NewSquareElement, indexes);
            SquareElements = [SquareElements; NewSquareElement]; %#ok 
            continue
        end
        % Store non-square elements obtained by cutting square elements
        % with boundary. Such non-square elements are not endowed with node
        % indexes yet.
        [NewElement, LocalSurfaceElements] = cut(NewSquareElement, fun, tol);
        NonSquareElements = [NonSquareElements; NewElement]; %#ok 
        SE = [SE; LocalSurfaceElements]; %#ok
    end
end

% AFTER ELIMINATING SQUARE ELEMENTS THAT ARE OUTSIDE DOMAIN, RE-DETERMINE
% INDEXES OF NODES USED BY SQUARE ELEMENTS
P = P(accepted_node,:);
acceptedindexes = zeros(Nrect,1);
acceptedindexes(accepted_node,1) = linspace(1,length(P),length(P))';
for i=1:length(SquareElements)
    setPind(SquareElements(i), acceptedindexes(SquareElements(i).Pind));
end

% DETERMINE SET OF NON-REPEATED NODES UP TO SMALL TOLERANCE
for i=1:length(NonSquareElements)
   P = [P; NonSquareElements(i).P]; %#ok 
end
P = uniquetol(P,tol,'ByRows',true);

% FIX ELEMENTS BY ASSIGNING NODE INDEXES AND ELIMINATING DUPLICATE NODES UP
% TO SMALL TOLERANCE
for i=1:length(SquareElements)
   [~, ind] = ismembertol(SquareElements(i).P,P,tol,'ByRows',true);
   setPind(SquareElements(i), ind);
   setP(SquareElements(i), P(ind,:));
end

for i=1:length(NonSquareElements)
   [~, ind] = ismembertol(NonSquareElements(i).P,P,tol,'ByRows',true);
   setPind(NonSquareElements(i), ind);
   setP(NonSquareElements(i), P(ind,:));
end

BulkElements = [SquareElements; NonSquareElements];

SurfElements = zeros(length(SE),2);
for i=1:length(SE)
    [~, ind] = ismembertol(SE{i},P,tol,'ByRows',true);
    SurfElements(i,:) = ind';
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

% CUTS GIVEN ELEMENT BY BOUNDARY OF DOMAIN
function [CutElement, LocalSurfaceElements] = cut(Element, fun, tol)
    P = Element.P;
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
        CutElement = [];
        return
    end
    indnewsort = sort(ind);
    Pnewsort = Pnew(indnewsort,:);
    CutElement = element2d(Pnewsort, false, false, [], mean(Pnewsort,1));
    indboundary = find(abs(fun(Pnewsort)) < tol);
    LocalSurfaceElements = cell(0,1);
    for i=1:length(indnewsort)-1
        if all(ismember(indnewsort([i i+1]), indboundary))
            LocalSurfaceElements = [LocalSurfaceElements; {Pnew(indnewsort([i i+1]),:)}]; %#ok
        end
    end
    if all(ismember(indnewsort([end 1]), indboundary))
            LocalSurfaceElements = [LocalSurfaceElements; {Pnew(indnewsort([end 1]),:)}];
    end
end

% COMPUTES INTERSECTION POINT BETWEEN AN EDGE AND THE BOUNDARY OF THE
% DOMAIN
function newP = intersectEdge(P1, P2, fun)
    fun_restricted = @(alpha) fun(alpha*P1 + (1-alpha)*P2);
    alpha = fzero(fun_restricted, [0 1]);
    newP = alpha*P1 + (1-alpha)*P2;
end