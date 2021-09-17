classdef element2ddummy
    %ELEMENT2DDUMMY is a bare-bone class that models a convex 2D element:
    % 1) contains minimal information on the element
    % 2) allows for element cutting in the unit sphere
    
    properties (SetAccess = private)
        P(:,3) % Vertexes
        issquare(1,1) logical
    end
    
    methods
        function obj = element2ddummy(P,issquare)     
            obj.P = P;
            obj.issquare = issquare;
        end
        
        function Earray = dummy2element(EDarray)
            for ei = 1:length(EDarray)
                if EDarray(ei).issquare
                    Earray(ei) = element2dsquare(EDarray(ei).P); %#ok
                else
                    Earray(ei) = element2d(EDarray(ei).P, mean(EDarray(ei).P)); %#ok
                end
            end
        end
        
        function CutElements = cut(obj)
            % Works only for cutting convex elements to fit in the sphere
            
            altered = false;
            
            for i=1:size(obj.P,1)
               if norm(obj.P(i,:)) > 1
                   altered = true;
                   obj.issquare = false;
                   obj.P(i,:) = obj.P(i,:)/norm(obj.P(i,:));
               end
            end
            
            if not(altered)
               CutElements = obj;
               return 
            end
            
            for i=1:size(obj.P,1)-2
                CutElements(i,1) = element2ddummy(obj.P([1 i+1 i+2],:), false); %#ok
            end
        end
    end
end

