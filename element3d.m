classdef element3d < element3dabstract
    % ELEMENT3D represents a VEM polyhedral element with k=1

    properties
        is_cube(1,1) logical
    end
    
    properties (SetAccess = private)
        
        % INITIALISED VIA CONSTRUCTOR
        Faces
        P % Vertexes
        Pind %Indexes of vertexes
        NVert
        NFaces
        P0(1,3) double % The element is star-shaped wrt P0
        
        % Computed by initElement
        Volume
        Centroid
        Diameter

        % Computed by getLocalMatrices
        K
        M
        C
    end
    
    properties (SetAccess = private, GetAccess = private)
        OutwardNormals(:,3) double
        hasMatrices = false;
    end
    
    methods (Access = private)
       function obj = initElement(obj)

           for i=1:obj.NFaces
                getLocalMatrices(obj.Faces(i));
           end

           if obj.is_cube
                % Compute centroid and P0
                obj.Centroid = mean(obj.P);
                obj.P0 = obj.Centroid;
           
                % Compute diameter
                obj.Diameter = 2*(norm(obj.P(1,:) - obj.Centroid));
           
                % Compute volume
                obj.Volume = obj.Faces(1).EdgeLength^3;
                
                return
           end
           
           % Compute volume and centroid
           volumes = zeros(obj.NFaces, 1);
           CV = [0 0 0];
           for i=1:length(obj.Faces)
               volumes(i) = abs(obj.Faces(i).OrientedArea * (obj.P0 - obj.Faces(i).P0)')/3;
               CV = CV + volumes(i)*(obj.P0 + 3*obj.Faces(i).Centroid)/4;
           end
           
           obj.Volume = sum(volumes);
           obj.Centroid = CV/obj.Volume;
           
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
       
       function obj = computeLocalMatrices(obj)
            % if obj.NVert == 4
            %     obj.M = (ones(4)+eye(4))*obj.Volume/20;
            %     obj.C = obj.M;
            %     B = [obj.P'; [1, 1, 1, 1]];
            %     Binv = inv(B);
            %     obj.K = Binv(:,1:3)*Binv(:,1:3)'*obj.Volume;
            %     return
            % end

            if obj.is_cube
                % DETERMINE ADJACENCY STRUCTURE
                adj_pairs = zeros(8);
                opp_face_pairs = zeros(8);
                opp_pairs = zeros(8);
            
                for i=1:8
                    for j=1:i-1
                        edgelen = obj.Faces(1).EdgeLength;
                        dist = norm(obj.P(i,:) - obj.P(j,:));
                        if abs(dist - edgelen) < edgelen/10
                            adj_pairs(i,j) = 1;
                            adj_pairs(j,i) = 1;
                        else
                            if abs(dist - obj.Faces(1).Diameter) < edgelen/10
                                opp_face_pairs(i,j) = 1;
                                opp_face_pairs(j,i) = 1;
                            else
                                opp_pairs(i,j) = 1;
                                opp_pairs(j,i) = 1;
                            end
                        end
                    end
                 end

                %COMPUTING LOCAL STIFFNESS MATRIX FROM B AND D (See Hitchhiker's)
                obj.K = obj.Diameter*(3*eye(8) + adj_pairs -   opp_face_pairs - 3*opp_pairs)*sqrt(3)/48 ...
                  + obj.Diameter*(2*eye(8) - adj_pairs + 0*opp_face_pairs +   opp_pairs)/4;

                %COMPUTING LOCAL MASS MATRIX FROM H,B AND D (See Hitchhiker's)
                obj.M = obj.Volume*(51*eye(8) - 22*adj_pairs + opp_face_pairs + 24*opp_pairs)/96;

                %COMPUTING LOCAL CONSISTENCY MATRIX FROM H,B AND D (See Hitchhiker's)
                obj.C = obj.Volume*(3*eye(8) + 2*adj_pairs + opp_face_pairs + 0*opp_pairs)/96;
            
                return
            end


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
            HP = H*PInablastar;
            PI0 = PInabla; %solo per k=1,2.
            obj.M = HP'*PInablastar + obj.Volume*(eye(obj.NVert)-PI0)'*(eye(obj.NVert)-PI0);
            obj.C = HP'*PInablastar;
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
        function obj = element3d(P, Faces, is_cube, Pind, P0)
            % ELEMENT3D Construct an instance of this class
            obj.P = P;
            obj.NVert = size(P,1);
            obj.Faces = Faces;
            obj.NFaces = length(Faces);
            if nargin >= 3
                obj.is_cube = is_cube;
            end
            if nargin >= 4
               obj.Pind = Pind; 
            end
            if nargin >= 5
               obj.P0 = P0;
            end
        end

        function newobj = copyElement3d(obj)
            if not(isequal(class(obj), 'element3d'))
                error('Wrong input class')
            end
            newobj = element3d(obj.P, obj.Faces, obj.is_cube, obj.Pind, obj.P0);
        end
        
        function P_ind_boundary = get_P_indexes_boundary(obj)
            P_ind_boundary = [];
            for i=1:length(obj.Faces)
                if obj.Faces(i).is_boundary
                   P_ind_boundary = [P_ind_boundary; obj.Faces(i).Pind]; %#ok
                end
            end
            P_ind_boundary = unique(P_ind_boundary);
        end

        function E = shiftElement(obj, v)
           newFaces = obj.Faces;
           for i=1:size(newFaces,1)
              newFaces(i) = shiftElement(obj.Faces(i), v); 
           end
           E =  element3d(obj.P+repmat(v,length(obj.P),1), newFaces, obj.is_cube,obj.Pind);
        end

        function EE = extrude(obj, Ncube)
            % Extrudes the faces that need to be extruded
            if not(obj.is_cube)
               error('Extrusion is only supported for cubic elements'); 
            end
            EE = [];
            I = eye(3);
            deltax = max(obj.P(:,1)) - min(obj.P(:,1));
            for i=1:6
                dir = obj.Faces(i).NormalDirection;
                fabsdir = abs(obj.Faces(i).P(1,dir));
                maxdir = max(obj.P(:,dir));
                mindir = min(obj.P(:,dir));
                maxabsdir = max(abs(obj.P(:,dir)));
                % discard faces that point toward the inside of the sphere
                % (3 out of 6)
                if maxdir*mindir > 0 && fabsdir < maxabsdir
                    continue
                end
                verse = sign(obj.Faces(i).P(1,dir) - mean(obj.P(:,dir)));
                % discard faces, based on normal direction, that are not to
                % be extruded
                if max(vecnorm((obj.Faces(i).P + I(dir,:)*deltax*verse)')) <=1
                   continue 
                end
                [~, indb] = ismember(obj.Faces(i).P,obj.P,'rows');
                obj.Faces(i).Pind = obj.Pind(indb);
                EE = [EE; extrude(obj.Faces(i), Ncube)]; %#ok
            end
        end

        function plot(obj, pert)
            if nargin == 1
               pert = 0; 
            end
           hold on
           for i=1:size(obj.Faces,1)
               PF = obj.Faces(i).P;
               % radii = vecnorm(PF')';
               fill3(PF(:,1), PF(:,2), PF(:,3), PF(:,1)*0 + 1 -pert - obj.is_cube);
           end
        end

        function plotSolution(obj, sol)
           hold on
           solE = sol(obj.Pind);
           for i=1:size(obj.Faces,1)
               PF = obj.Faces(i).P;
               [tf,loc] = ismember(obj.P,PF,'rows');
               [~,p] = sort(loc(tf));
               indF = find(tf);
               indF = indF(p);
               if not(isempty(indF))
                    fill3(PF(:,1), PF(:,2), PF(:,3), solE(indF),'FaceColor', 'interp', 'EdgeColor', 'none');
               end
           end
        end

        function obj = getLocalMatrices(obj)
            if not(obj.hasMatrices)
               obj = initElement(obj);
               obj = computeLocalMatrices(obj);
               obj.hasMatrices = true;
            end
        end

        function setP(obj, P)
            obj.P = P;
        end

        function setPind(obj, Pind)
            obj.Pind = Pind;
        end
    end
end

