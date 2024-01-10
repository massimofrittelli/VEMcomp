classdef (Abstract) element2dabstract < handle
    %element2dabstract is an abstract class that defines the functionality of
    %any 2D VEM element
    
    properties (Abstract, SetAccess = private)
        P(:,3) double % Nodes
        P0(1,3) double % Element is star-shaped wrt P0
        NVert(1,1) double % Numver of nodes
        Area(1,1) double
        OrientedArea(1,3) double
        Centroid (1,3) double
        Diameter(1,1) double
        K(:,:) double % Stiffness matrix
        M(:,:) double % Mass matrix
        C(:,:) double % Consistency matrix
    end
    
            
end

