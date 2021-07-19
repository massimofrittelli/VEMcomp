classdef element2dsquare < element2d
    % ELEMENT2D represents a square VEM element with k=1
    
    properties
        EdgeLength(1,1) double
    end
    
    methods
        function obj = element2dsquare(edgeLength)
            % ELEMENT2DSQUARE Construct an instance of this class
            obj@element2d([]);
            obj.EdgeLength = edgeLength;
        end
        
        function [KE, ME] = local_matrices(obj)
            % LOCALMATRICES computes local stiffness and mass matrices
            % Dummy implementation
            KE = diag(obj.Nodes(:,1));
            ME = diag(obj.Nodes(:,2));
        end
    end
end

