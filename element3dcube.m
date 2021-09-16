classdef element3dcube < element3dabstract
    % ELEMENT3DCUBE represents a cubic VEM element with k=1
        
    properties (SetAccess = private)
        
        % INITIALISED VIA CONSTRUCTOR
        Faces
        P % Vertices
        
        % COMPUTED UPON INITIALISATION
        NVert
        NFaces
        Volume
        Centroid
        Diameter
        K
        M
    end
    
    methods (Access = private)
       function obj = initElement(obj)
           
           % Compute centroid
           obj.Centroid = mean(obj.P);
           
           % Compute number of vertices and faces
           obj.NVert = 8;
           obj.NFaces = 6;
           
           % Compute diameter
           obj.Diameter = 2*(norm(obj.P(1,:) - obj.Centroid));
           
           % Compute volume
           obj.Volume = obj.Faces(1).EdgeLength^3;
           
       end
       
       function obj = setLocalMatrices(obj)
           
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
            
       end
       
    end
    
    methods
        function obj = element3dcube(Faces, P)
            % ELEMENT3D Construct an instance of this class
            obj.Faces = Faces;
            obj.P = P;
            obj = initElement(obj);
            obj = setLocalMatrices(obj);
        end
        
    end
end

