classdef element2d < element2dabstract
    % ELEMENT2D represents a VEM polygonal element with k=1
    
    properties(SetAccess = private)
        
        % CONSTRUCTOR INPUTS
        P % Vertexes
        P0 % Element is star-shaped wrt P0
        Pind % Indexes of vertexes
        is_boundary
        
        % COMPUTED BY CONSTRUCTOR
        NVert
        Area
        OrientedArea
        Centroid
        Diameter
        K
        M
        CM % Consistency matrix
    end
    
    properties(SetAccess = private, GetAccess = private)
        TransformedP(:,2) double
        TransformedP0(1,2) double
        TransformedCentroid(1,2) double
    end
    
    methods(Access = private)
        function obj = initElement(obj)
            
             % Compute number of vertices
             obj.NVert = length(obj.P);
             
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
        end
        
        function obj = setLocalMatrices(obj)
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
            C = H*PInablastar;
            PI0 = PInabla; %solo per k=1,2.
            obj.CM = C'*PInablastar;
            obj.M = C'*PInablastar + obj.Area*(eye(obj.NVert)-PI0)'*(eye(obj.NVert)-PI0);
        end
        
        
        function [XY,W] = quadrature_quadratic(obj)
            % Determines quadrature weights and nodes on polygon
            % W = zeros(4*obj.NVert,1);
            % XY = zeros(4*obj.NVert,2);
            % for i = 1:obj.NVert-1
            %     PP = [obj.TransformedP(i:i+1,:); obj.TransformedP0];
            %     [XY(4*i-3:4*i,:), W(4*i-3:4*i,1)] = quadrature_triangle_quadratic(PP);
            % end
            % i = obj.NVert;
            % PP = [obj.TransformedP([obj.NVert,1],:); obj.TransformedP0];
            % [XY(4*i-3:4*i,:), W(4*i-3:4*i,1)] = quadrature_triangle_quadratic(PP);
            
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
        function obj = element2d(P, P0, Pind, is_boundary)
            % ELEMENT Construct an instance of this class
            obj.P = P;
            obj.P0 = P0;
            if nargin >= 3
               obj.Pind = Pind;
            end
            if nargin >= 4
               obj.is_boundary = is_boundary;
            end
            obj = initElement(obj);
            obj = setLocalMatrices(obj);
        end
    end
end

