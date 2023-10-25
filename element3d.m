classdef element3d < element3dabstract
    % ELEMENT3D represents a VEM polyhedral element with k=1
    
    properties (SetAccess = private)
        
        % INITIALISED VIA CONSTRUCTOR
        Faces
        P % Vertexes
        Pind %Indexes of vertexes
        
        % INITIALISED VIA CONSTRUCTOR - SUBCLASS SPECIFIC
        P0(1,3) double % The element is star-shaped wrt P0
        
        % COMPUTED UPON INITIALISATION
        NVert
        NFaces
        Volume
        Centroid
        Diameter
        K
        M
        CM
    end
    
    properties (SetAccess = private, GetAccess = private)
        OutwardNormals(:,3) double
    end
    
    methods (Access = private)
       function obj = initElement(obj)
           
           % Compute number of vertices and faces
           obj.NVert = length(obj.P);
           obj.NFaces = length(obj.Faces);
           
           % Compute volume and centroid
           volumes = zeros(obj.NFaces, 1);
           C = [0 0 0];
           for i=1:length(obj.Faces)
               volumes(i) = abs(obj.Faces(i).OrientedArea * (obj.P0 - obj.Faces(i).P0)')/3;
               C = C + volumes(i)*(obj.P0 + 3*obj.Faces(i).Centroid)/4;
           end
           
           obj.Volume = sum(volumes);
           obj.Centroid = C/obj.Volume;
           
           % Compute diameter
           distances = zeros(obj.NVert);
           for i=1:obj.NVert
              for j=i+1:obj.NVert
                  distances(i,j) = norm(obj.P(i,:)-obj.P(j,:)); 
              end
           end 
             
           obj.Diameter = max(max(distances));
           
           % Compute OUTWARD normals on each face
           ON = zeros(obj.NFaces,3);
           for i=1:obj.NFaces
              ON(i,:) = obj.Faces(i).OrientedArea/norm(obj.Faces(i).OrientedArea);
              ON(i,:) = ON(i,:)*sign((obj.Faces(i).P0 - obj.P0)*ON(i,:)');
           end
           
           obj.OutwardNormals = ON;
           
       end
       
       function obj = setLocalMatrices(obj)
           % computing gradient of the non-constant barycentric monomials 
            nabla2 = [1,0,0]/obj.Diameter;
            nabla3 = [0,1,0]/obj.Diameter;
            nabla4 = [0,0,1]/obj.Diameter;
            
            % computing outward normal derivatives of the non-constant
            % barycentric monomials on the faces
            
            normders2 = obj.OutwardNormals * nabla2';
            normders3 = obj.OutwardNormals * nabla3';
            normders4 = obj.OutwardNormals * nabla4';

            %computing matrix B (see Hitchhiker's)
            B = zeros(4,obj.NVert);
            B(1,:) = ones(1,obj.NVert)/obj.NVert;
            for i=1:obj.NVert
                B(2,i) = boundaryIntegral(obj, normders2, i);
                B(3,i) = boundaryIntegral(obj, normders3, i);
                B(4,i) = boundaryIntegral(obj, normders4, i);
            end

            %barycentric monomials
            monomials = {@(x,y,z) 0*x(:,1) + 1;
                @(x) (x(:,1)-obj.Centroid(1))/obj.Diameter;
                @(x) (x(:,2)-obj.Centroid(2))/obj.Diameter;
                @(x) (x(:,3)-obj.Centroid(3))/obj.Diameter};
         
            %computing matrix D (see Hitchhiker's)
            D = zeros(obj.NVert,4);
            D(:,1) = ones(obj.NVert,1);
            for j=2:4
                for i=1:obj.NVert
                    D(i,j) = monomials{j}(obj.P(i,:));
                end
            end

            %COMPUTING LOCAL STIFFNESS MATRIX FROM B AND D (See Hitchhiker's)
            G = B*D;
            Gtilde = G;
            Gtilde(1,:) = zeros(1,4);
            PInablastar = G\B;
            PInabla = D*PInablastar;
            obj.K = PInablastar'*Gtilde*PInablastar + obj.Diameter*(eye(obj.NVert)-PInabla)'*(eye(obj.NVert)-PInabla);

            %computing matrix H (see Hitchhiker's)
            [XYZ,W] = quadrature(obj);
            H = zeros(4,4);
            for i=1:4
                for j=1:4
                    fun = @(x) monomials{i}(x).*monomials{j}(x);
                    H(i,j) = W'*fun(XYZ);
                end
            end

            %COMPUTING LOCAL MASS MATRIX FROM H,B AND D (See Hitchhiker's)
            C = H*PInablastar;
            PI0 = PInabla; %solo per k=1,2.
            obj.M = C'*PInablastar + obj.Volume*(eye(obj.NVert)-PI0)'*(eye(obj.NVert)-PI0);
            obj.CM = C'*PInablastar;
       end
       
       function [XYZ,W] = quadrature(obj)
            % COMPUTES QUADRATURE NODES AND WEIGHTS ON THE WHOLE 3D ELEMENT
            XYZ = []; W = []; 
            for i = 1:obj.NFaces
               [XYZF, WF] = quadraturePyramid(obj, obj.Faces(i));
               XYZ = [XYZ; XYZF]; W = [W;WF];  %#ok
            end
       end
        
       function [XYZ,W] = quadraturePyramid(obj, face)
           % COMPUTES QUADRATURE NODES AND WEIGHTS ON PYRAMID HAVING AN
           % element2d AS BASE
           XYZ = zeros(4*face.NVert,3); W = zeros(4*face.NVert,1); 
           for i=1:face.NVert-1
              PP = [face.P([i i+1],:); face.P0; obj.P0];
              [XYZ(4*i-3:4*i,:), W(4*i-3:4*i,:)] = quadrature_tetrahedron_quadratic(PP);
           end
           i=face.NVert;
           PP = [face.P([face.NVert 1],:); face.P0; obj.P0];
           [XYZ(4*i-3:4*i,:), W(4*i-3:4*i,:)] = quadrature_tetrahedron_quadratic(PP);
       end
       
       function I = boundaryIntegral(obj, normders, i)
           I = 0;
           for j=1:obj.NFaces
              ii = find(ismember(obj.Faces(j).P, obj.P(i,:),'rows'), 1);
              if not(isempty(ii))
                  I = I + normders(j)*sum(obj.Faces(j).M(ii,:));
              end
           end
       end
    end
    
    methods
        function obj = element3d(Faces, P, P0, Pind)
            % ELEMENT3D Construct an instance of this class
            obj.Faces = Faces;
            obj.P = P;
            obj.P0 = P0;
            if nargin >= 4
               obj.Pind = Pind; 
            end
            obj = initElement(obj);
            obj = setLocalMatrices(obj);
        end
        
    end
end

