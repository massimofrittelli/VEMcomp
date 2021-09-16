classdef element2dsquare < element2dabstract
    % ELEMENT2DSQUARE represents a VEM square element with k=1
    
    properties(SetAccess = private)
        
        % CONSTRUCTOR INPUT
        P % Nodes
        
        % COMPUTED BY CONSTRUCTOR
        P0 % Element is star-shaped wrt P0
        NVert
        Area
        OrientedArea
        Centroid
        Diameter
        K
        M
        
        % COMPUTED BY CONSTRUCTOR AND NOT PRESENT IN ABSTRACT CLASS
        EdgeLength(1,1) double
    end
    
    methods(Access = private)
        
        function obj = initElement(obj)
            
             % Compute number of vertices
             obj.NVert = 4;
             
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
             
        end
        
        function obj = setLocalMatrices(obj)
            
            % Compute local stiffness matrix K
            obj.K = (4*eye(4) - ones(4,4))/4;
            
            % Compute local mass matrix M
            v2 = ones(2,1);
            v3 = ones(3,1);
            obj.M = obj.Area*(17*eye(4) -9*(diag(v3,1) + diag(v3,-1)) + 13*(diag(v2,2) + diag(v2,-2)) -9*(diag(1,3) + diag(1,-3)))/48;
        end
        
    end
    
    methods
        function obj = element2dsquare(P)
           % ELEMENT2DSQUARE Construct an instance of this class
            obj.P = P;
            obj = initElement(obj);
            obj = setLocalMatrices(obj);
        end
    end
end

