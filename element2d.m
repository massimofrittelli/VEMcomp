classdef element2d < handle
    % ELEMENT2D represents a VEM polygonal element with k=1

    properties
        to_plot logical
    end
    
    properties(SetAccess = private)
        
        %Initialized by constructor
        P(:,3) double % Nodes % Vertexes
        P0(1,3) double % Element is star-shaped wrt P0
        Pind % Indexes of vertexes
        is_square
        is_boundary
        NVert(1,1) double
        
        %Computed by getLocalMatrices
        Area(1,1) double
        OrientedArea(1,3) double
        Centroid(1,3) double
        Diameter(1,1) double
        EdgeLength(1,1) double
        K(:,:) double % Stiffness matrix
        M(:,:) double % Mass matrix
        C(:,:) double % Consistency matrix       
    end
    
    properties(SetAccess = private, GetAccess = private)
        TransformedP(:,2) double
        TransformedP0(1,2) double
        TransformedCentroid(1,2) double
        hasMatrices = false;
    end

    properties (Dependent)
       NormalDirection(3,1) double 
    end
    
    methods(Access = private)
        
        function getLocalMatricesSquare(obj)
                % Compute edge length
                obj.EdgeLength = norm(obj.P(1,:) - obj.P(2,:));
             
                % Compute Area
                obj.Area = obj.EdgeLength^2;
             
                % Compute Oriented Area
                obj.OrientedArea = cross(obj.P(3,:) - obj.P(2,:), obj.P(2,:) - obj.P(1,:));
             
                % Compute Centroid and P0
                obj.Centroid = (obj.P(1,:) + obj.P(3,:))/2;
                obj.P0 = obj.Centroid;
             
                % Compute Diameter
                obj.Diameter = norm(obj.P(1,:) - obj.P(3,:));

                % Compute local stiffness matrix K
                obj.K = (4*eye(4) - ones(4,4))/4;
            
                % Compute local mass matrix M
                v2 = ones(2,1);
                v3 = ones(3,1);
                obj.M = obj.Area*(17*eye(4) -9*(diag(v3,1) + diag(v3,-1)) + 13*(diag(v2,2) + diag(v2,-2)) -9*(diag(1,3) + diag(1,-3)))/48;
        
                % Compute consistency matrix CM
                obj.C = obj.Area*(5*eye(4) +3*(diag(v3,1) + diag(v3,-1)) + 1*(diag(v2,2) + diag(v2,-2)) +3*(diag(1,3) + diag(1,-3)))/48;
            
        end
        
        function getLocalMatricesPoly(obj)

            % if obj.NVert == 3
            %     obj.M = (ones(3)+eye(3))*obj.Area/12;
            %     obj.C = obj.M;
            %     % edges of the considered triangle
            %     v1 = obj.P(3,:)-obj.P(2,:);
            %     v2 = obj.P(1,:)-obj.P(3,:);
            %     v3 = obj.P(2,:)-obj.P(1,:);
            %     % heights of the considered triangle
            %     h = [v2-(v1*v2')*v1/norm(v1)^2;
            %          v3-(v2*v3')*v2/norm(v2)^2;
            %          v1-(v3*v1')*v3/norm(v3)^2];
            %     % gradients of the local basis functions
            %     nabla = [h(1,:)/norm(h(1,:))^2;
            %              h(2,:)/norm(h(2,:))^2;
            %              h(3,:)/norm(h(3,:))^2];
            %     obj.K = nabla*nabla'*obj.Area*2;
            %     return
            % end

            % Compute area, oriented area and centroid
             orientedAreas = zeros(obj.NVert, 3);
             areas = zeros(obj.NVert, 1);
             centroids = zeros(obj.NVert,3);
             for i=1:obj.NVert-1
                orientedAreas(i,:) = cross((obj.P(i,:)-obj.P0), (obj.P(i+1,:)-obj.P0));
                areas(i) = norm(orientedAreas(i,:));
                centroids(i,:) = (obj.P(i,:) + obj.P(i+1,:) + obj.P0);
             end
             orientedAreas(obj.NVert,:) = cross((obj.P(obj.NVert,:)-obj.P0), (obj.P(1,:)-obj.P0));
             areas(obj.NVert) = norm(orientedAreas(obj.NVert,:));
             centroids(obj.NVert,:) = (obj.P(obj.NVert,:) + obj.P(1,:) + obj.P0);
             
             obj.OrientedArea = sum(orientedAreas,1)/2;
             obj.Area = sum(areas)/2;
             obj.Centroid = sum(centroids.*repmat(areas,[1, 3]),1)/(3*sum(areas));
             
             % Compute diameter
             distances = zeros(obj.NVert);
             for i=1:obj.NVert
                for j=i+1:obj.NVert
                    distances(i,j) = norm(obj.P(i,:)-obj.P(j,:)); 
                end
             end 
             
             obj.Diameter = max(max(distances));
             
             % Compute element transformed to xy-plane
             % 1) shifting element by moving first node PE(1,:) to the origin
             PS = obj.P - repmat(obj.P(1,:),obj.NVert,1); % shifted nodes
             GS = obj.Centroid - obj.P(1,:); % shifted centroid
             P0S = obj.P0 - obj.P(1,:); % shifted star-shaped point

             % 2) rotating shifted element to the (x,y) plane
             firstedge = PS(2,:);
             normfirstedge = norm(firstedge);
             xPR = zeros(obj.NVert,1); % x's of rotated nodes
             yPR = zeros(obj.NVert,1); % y's of rotated nodes
             GR = zeros(1,2); % rotated centroid
             P0R = zeros(1,2); % rotated star-shaped point
             for i=2:obj.NVert
                 xPR(i) = PS(i,:)*firstedge'/normfirstedge;
                 yPR(i) = norm(cross(PS(i,:),firstedge))/normfirstedge;
             end
             PR = [xPR,yPR];
             GR(1) = GS*firstedge'/normfirstedge;
             GR(2) = norm(cross(GS,firstedge))/normfirstedge;
             P0R(1) = P0S*firstedge'/normfirstedge;
             P0R(2) = norm(cross(P0S,firstedge))/normfirstedge;
             
             obj.TransformedP = PR;
             obj.TransformedP0 = P0R;
             obj.TransformedCentroid = GR;

            %computing edges of the element
            edges = zeros(obj.NVert,2);
            for i=1:obj.NVert-1
                edges(i,:) = obj.TransformedP(i+1,:)-obj.TransformedP(i,:); 
            end
            edges(obj.NVert,:) = obj.TransformedP(1,:)-obj.TransformedP(end,:);

            %computing gradient of the non-constant scaled monomials 
            nabla1 = [1,0]/obj.Diameter;
            nabla2 = [0,1]/obj.Diameter;

            %computing unit (outward ?) normals to each edge
            normals = -edges*[0,-1;1,0]';
            for i=1:obj.NVert
                normals(i,:) = normals(i,:)/norm(normals(i,:)); 
            end

            %computing outward normal derivatives of the non-constant scaled monomials
            %along the edges
            normders1 = normals * nabla1';
            normders2 = normals * nabla2';

            %computing matrix B (see Hitchhiker's)
            B = zeros(3,obj.NVert);
            B(1,:) = ones(1,obj.NVert)/obj.NVert;
            B(2,1) = (normders1(obj.NVert)*norm(edges(obj.NVert,:))+normders1(1)*norm(edges(1,:)))/2;
            B(3,1) = (normders2(obj.NVert)*norm(edges(obj.NVert,:))+normders2(1)*norm(edges(1,:)))/2;
            for i=2:obj.NVert
                B(2,i) = (normders1(i-1)*norm(edges(i-1,:))+normders1(i)*norm(edges(i,:)))/2;
                B(3,i) = (normders2(i-1)*norm(edges(i-1,:))+normders2(i)*norm(edges(i,:)))/2;
            end

            %barycentric monomials
            monomials = {@(x,y) 0*x + 1;
                @(x,y) (x-obj.TransformedCentroid(1))/obj.Diameter;
                @(x,y) (y-obj.TransformedCentroid(2))/obj.Diameter};
         
            %computing matrix D (see Hitchhiker's)
            D = zeros(obj.NVert,3);
            D(:,1) = ones(obj.NVert,1);
            for j=2:3
                for i=1:obj.NVert
                    D(i,j) = monomials{j}(obj.TransformedP(i,1),obj.TransformedP(i,2));
                end
            end

            %COMPUTING LOCAL STIFFNESS MATRIX FROM B AND D (See Hitchhiker's)
            G = B*D;
            Gtilde = G;
            Gtilde(1,:) = zeros(1,3);
            PInablastar = G\B;
            PInabla = D*PInablastar;
            obj.K = PInablastar'*Gtilde*PInablastar + (eye(obj.NVert)-PInabla)'*(eye(obj.NVert)-PInabla);

            %computing matrix H (see Hitchhiker's)
            [XY,W] = quadrature_quadratic(obj);
            H = zeros(3,3);
            for i=1:3
                for j=1:3
                    fun = @(x,y) monomials{i}(x,y).*monomials{j}(x,y);
                    H(i,j) = W' * fun(XY(:,1),XY(:,2));
                end
            end

            %COMPUTING LOCAL MASS MATRIX FROM H,B AND D (See Hitchhiker's)
            HP = H*PInablastar;
            PI0 = PInabla; %solo per k=1,2.
            obj.C = HP'*PInablastar;
            obj.M = HP'*PInablastar + obj.Area*(eye(obj.NVert)-PI0)'*(eye(obj.NVert)-PI0);
        end
        
        
        function [XY,W] = quadrature_quadratic(obj)
            W = zeros(3*obj.NVert,1);
            XY = zeros(3*obj.NVert,2);
            for i = 1:obj.NVert-1
                PP = [obj.TransformedP(i:i+1,:); obj.TransformedP0];
                [XY(3*i-2:3*i,:), W(3*i-2:3*i,1)] = quadrature_triangle_quadratic(PP);
            end
            i = obj.NVert;
            PP = [obj.TransformedP([obj.NVert,1],:); obj.TransformedP0];
            [XY(3*i-2:3*i,:), W(3*i-2:3*i,1)] = quadrature_triangle_quadratic(PP);
        end
        
    end
    
    methods
        function obj = element2d(P, is_square, is_boundary, Pind, P0)
            % ELEMENT Construct an instance of this class
            obj.P = P;
            obj.NVert = length(obj.P);
            if nargin >= 2
                obj.is_square = is_square;
            end
            if nargin >= 3
               obj.is_boundary = is_boundary;
            end
            if nargin >= 4
               obj.Pind = Pind;
            end
            if nargin >= 5
                obj.P0 = P0;
            end
        end

        function newobj = copyElement2d(obj)
            if not(isequal(class(obj), 'element2d'))
                error('Wrong input class')
            end
            newobj = element2d(obj.P, obj.is_square, obj.is_boundary, obj.Pind, obj.P0);
        end


        function ND = get.NormalDirection(obj)
           if not(obj.is_square)
              error('Normal direction is supported square elements only') 
           end
           direction = cross(obj.P(3,:)-obj.P(2,:) , obj.P(2,:)-obj.P(1,:));
           [~, ND] = max(abs(direction));
        end

        function E = shiftElement(obj, v)
           E =  element2d(obj.P+repmat(v,length(obj.P),1),obj.is_square,false,obj.Pind);
        end

        function EE = extrude(obj, Ncube)
           EP = zeros(size(obj.P));
           actuallyExtruded = logical(size(obj.P,1));
           for i=1:obj.NVert
               [EP(i,:), actuallyExtruded(i)] = extrude_node(obj.P(i,:));
           end
           extruded_ind = obj.Pind + Ncube;
           extruded_ind(not(actuallyExtruded)) = obj.Pind(not(actuallyExtruded));
           ExtrudedFaces = obj;
           
           % Re-ordering vertexes of square face of extruded element, in
           % such a wat that the normal is outward
           ExtrudedFaces.P = flipud(ExtrudedFaces.P);
           ExtrudedFaces.Pind = flipud(ExtrudedFaces.Pind);
           
           % CREATE TRIANGULAR FACES LYING ON THE SURFACE
           if norm(EP(1,:) + EP(3,:)) >= norm(EP(2,:) + EP(4,:))
               NewExtrudedFace1 = element2d(EP([1 2 3],:), false, true);
               NewExtrudedFace1.Pind = extruded_ind([1 2 3],1);
               NewExtrudedFace2 = element2d(EP([1 3 4],:), false, true);
               NewExtrudedFace2.Pind = extruded_ind([1 3 4],1);
           else
               NewExtrudedFace1 = element2d(EP([1 2 4],:), false, true);
               NewExtrudedFace1.Pind = extruded_ind([1 2 4],1);
               NewExtrudedFace2 = element2d(EP([2 3 4],:), false, true);
               NewExtrudedFace2.Pind = extruded_ind([2 3 4],1);
           end
               ExtrudedFaces = [ExtrudedFaces; NewExtrudedFace1; NewExtrudedFace2];
           % CREATE LATERAL SQUARE FACES
           for i=1:obj.NVert
              LP = unique([obj.P([i 1+rem(i,obj.NVert)],:); EP([1+rem(i,obj.NVert) i],:)],'rows','stable');
              NewExtrudedFace = element2d(LP, false, false);
              EPind = unique([obj.Pind([i 1+rem(i,obj.NVert)],1); extruded_ind([1+rem(i,obj.NVert) i],1)],'stable');
              NewExtrudedFace.Pind = EPind;
              ExtrudedFaces = [ExtrudedFaces; NewExtrudedFace]; %#ok 
           end
           EE = element3d([obj.P; EP(actuallyExtruded,:)], ExtrudedFaces, false, [obj.Pind; extruded_ind(actuallyExtruded)]);
        end

        function plot(obj, faceColor, edgeColor, edgeAlpha)
            switch nargin
                case 1
                    fill3(obj.P(:,1), obj.P(:,2), obj.P(:,3), obj.P(:,1)*0 + 1 - obj.is_square);
                case 2
                    fill3(obj.P(:,1), obj.P(:,2), obj.P(:,3), faceColor);
                case 3
                    fill3(obj.P(:,1), obj.P(:,2), obj.P(:,3), faceColor, 'EdgeColor', edgeColor);
                case 4
                    fill3(obj.P(:,1), obj.P(:,2), obj.P(:,3), faceColor, 'EdgeColor', edgeColor, 'EdgeAlpha', edgeAlpha);
            end
        end

        function obj = getLocalMatrices(obj)
            if obj.hasMatrices
                return
            end

            if obj.is_square
                getLocalMatricesSquare(obj);
                obj.hasMatrices = true;
                return;
            end

            getLocalMatricesPoly(obj);
            obj.hasMatrices = true;    
        end

        function setP(obj, P)
            obj.P = P;
        end

        function setPind(obj, Pind)
            obj.Pind = Pind;
        end
    end
end

function [P, wasActuallyExtruded] = extrude_node(P)
    tol = 1e-8;
    % TYPE 1: LESS UNIFORM, BUT INTERNAL FACES ARE GUARANTEED TO BE FLAT
    n = signtol(P);
    alpha = (- P*n' + sqrt((P*n')^2 -norm(n)^2*(norm(P)^2-1)))/norm(n)^2;
    wasActuallyExtruded = (abs(alpha) > tol);
    P = P + wasActuallyExtruded*alpha*n;
    % TYPE 2: MORE UNIFORM, BUT INTERNAL FACES ARE NOT GUARANTEED TO BE
    % FLAT
%     wasActuallyExtruded = (1 - norm(P) > tol);
%     P = P/norm(P);
end

function x = signtol(x)
    tol = 1e-10;
    x = sign(x).*(abs(x) > tol);
end

