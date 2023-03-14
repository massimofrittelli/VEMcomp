function [P, h, SquareElements, NonSquareElements] = generate_mesh_flat_domain(fun, xmax, Nx)
%GENERATE_MESH_FLAT_DOMAIN Summary of this function goes here
%   Detailed explanation goes here

Nsquare = Nx^2;
hx = 2*xmax/(Nx-1); % Discretisation step along each dimension
h = hx*sqrt(3); % Meshsize

% GENERATE ONE SQUARE ELEMENT
PS = [0 0 0; 0 1 0; 1 1 0; 1 0 0]*hx;
ES = element2d_dummy(PS, true);

% CREATING GRIDPOINTS OF BOUNDING BOX
x = linspace(-xmax,xmax,Nx); % Gridpoints in [-1,1]
P = zeros(Nsquare,3); % Gridpoints of bounding box
for i=0:Nx-1
    for j=0:Nx-1  
    	P(Nx*i+j+1,:) = [x(i+1) x(j+1) 0];
    end
end

% GENERATE ELEMENTS
SquareElements = [];
NonSquareElements = [];
accepted_node = false(Nsquare,1);
for i=0:Nx-2 % For each element of the bounding box
    for j=0:Nx-2
        indexes = [Nx*i+j+1
                   Nx*i+j+2
                   Nx*(i+1)+j+2
                   Nx*(i+1)+j+1];                  
        NewSquareElement = shiftElement(ES, P(indexes(1),:));
        if is_outside(NewSquareElement, fun)
            continue
        end
        if is_inside(NewSquareElement, fun)
            % Store square elements that are inside domain
            accepted_node(indexes) = true(4,1);
            NewSquareElement.Pind = indexes;
            SquareElements = [SquareElements; NewSquareElement]; %#ok 
            continue
        end
        % Store non-square elements obtained by cutting square elements
        % with boundary. Such non-square elements are not endowed with node
        % indexes yet.
        NewElement = cut(NewSquareElement, fun);
        NonSquareElements = [NonSquareElements; NewElement]; %#ok 
    end
end

% AFTER ELIMINATING SQUARE ELEMENTS THAT ARE OUTSIDE DOMAIN, RE-DETERMINE
% INDEXES OF NODES USED BY SQUARE ELEMENTS
P = P(accepted_node,:);
acceptedindexes = zeros(Nsquare,1);
acceptedindexes(accepted_node,1) = linspace(1,length(P),length(P))';
for i=1:length(SquareElements)
   SquareElements(i).Pind = acceptedindexes(SquareElements(i).Pind); %#ok
end

for i=1:length(NonSquareElements)
   P = [P; NonSquareElements(i).P]; %#ok 
end
P = unique(P,'rows','stable');

for i=1:length(NonSquareElements)
   [~,ind1, ind2] = intersect(P,NonSquareElements(i).P,'rows','legacy');
   NonSquareElements(i).Pind = nonzeros(ind1(ind2)); %#ok 
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
function CutElement = cut(Element, fun)
    P = Element.P;
    strictly_inside = fun(P) < 0;
    Pnew = [];
    for i=1:length(P)
       j = rem(i, length(P)) + 1;
       if strictly_inside(i)
           Pnew = [Pnew; P(i,:)]; %#ok
       end
       if sum(strictly_inside([i j])) == 1
           Pnew = [Pnew; intersectEdge(P(i,:), P(j,:), fun)]; %#ok
       end
    end
    CutElement = element2d_dummy(Pnew, false);
end

% COMPUTES INTERSECTION POINT BETWEEN AN EDGE AND THE BOUNDARY OF THE
% DOMAIN
function newP = intersectEdge(P1, P2, fun)
    fun_restricted = @(alpha) fun(alpha*P1 + (1-alpha)*P2);
    alpha = fzero(fun_restricted, [0 1]);
    newP = alpha*P1 + (1-alpha)*P2;
end