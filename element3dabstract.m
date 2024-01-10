classdef (Abstract) element3dabstract < handle
    %element3dabstract is an abstract class that defines the functionality of
    %any 3D VEM element
    
    properties (Abstract, SetAccess = private)
        
        Faces(:,1) element2d
        P(:,3) double % Vertices
        NVert(1,1) double % Number of vertices
        NFaces(1,1) double % Number of faces
        Volume(1,1) double
        Centroid(1,3) double
        Diameter(1,1) double
        K(:,:) double % Stiffness matrix
        M(:,:) double % Mass matrix
    end
    
            
end

