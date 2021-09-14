classdef (Abstract) element3dgeneric < matlab.mixin.Heterogeneous
    %ELEMENT3DGENERIC is an absract class that defines the functionality of
    %any 3D VEM element
    
    properties (Abstract, SetAccess = private)
        
        Faces(:,1) element2dgeneric
        P(:,3) double % Vertices
        P0(1,3) double % The element is star-shaped wrt P0
        NVert(1,1) double % Number of vertices
        NFaces(1,1) double % Number of faces
        Volume(1,1) double
        Centroid(1,3) double
        Diameter(1,1) double
        K(:,:) double % Stiffness matrix
        M(:,:) double % Mass matrix
    end
    
            
end

