classdef element3dcube < element3d
    % ELEMENT3DCUBE represents a cubic VEM element with k=1
    
    properties
        EdgeLength(1,1) double
    end
    
    methods
        function obj = element3dcube(edgeLength)
            % ELEMENT3DCUBE Construct an instance of this class
            obj@element3d(8,6,[]);
            obj.EdgeLength = edgeLength;
        end
        
        function [KE, ME] = localMatrices(obj)
            % LOCALMATRICES computes local stiffness and mass matrices
            % Dummy implementation
            KE = obj.Diameter*eye(8);
            ME = obj.Diameter*eye(8);
        end
    end
end

