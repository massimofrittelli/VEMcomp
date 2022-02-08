classdef element2ddummy
    %ELEMENT2DDUMMY is a bare-bone class that models a convex 2D element:
    % 1) contains minimal information on the element
    % 2) allows for element extrusion in the unit sphere
    % 3) allows for element translation/shifting
    
    properties
        Pind(:,1) double = [] % Indexes of vertexes
    end
    
    properties (SetAccess = private)
        P(:,3) double % Vertexes
        NVert(1,1) double % Number of vertexes
        issquare(1,1) logical
    end
    
    properties (Dependent)
       NormalDirection(3,1) double 
    end
    
    methods
        function obj = element2ddummy(P,issquare)     
            obj.P = P;
            obj.NVert = size(P,1);
            obj.issquare = issquare;
        end
        
        function ND = get.NormalDirection(obj)
           if not(obj.issquare)
              error('Normal direction is supported square elements only') 
           end
           direction = cross(obj.P(3,:)-obj.P(2,:) , obj.P(2,:)-obj.P(1,:));
           [~, ND] = max(abs(direction));
        end
        
        function E = shiftElement(obj, v)
           E =  element2ddummy(obj.P+repmat(v,length(obj.P),1),obj.issquare);
        end
        
        function EE = extrude(obj,Pind, Ncube)
           EP = obj.P./repmat(vecnorm(obj.P')',[1 3]); % Extruded points
           ExtrudedFaces = obj;
           for i=1:obj.NVert-2
               NewExtrudedFace = element2ddummy(EP([1 i+1 i+2],:), false);
               if not(isempty(obj.Pind))
                    NewExtrudedFace.Pind = obj.Pind([1 i+1 i+2],1) + Ncube;
               end
               ExtrudedFaces = [ExtrudedFaces; NewExtrudedFace]; %#ok
           end
           for i=1:obj.NVert
              EP1 = [obj.P([i 1+rem(i,obj.NVert)],:); EP(i,:)];
              EP2 = [EP([i 1+rem(i,obj.NVert)],:); obj.P(1+rem(i,obj.NVert),:)];
              NewExtrudedFace1 = element2ddummy(EP1, false);
              NewExtrudedFace2 = element2ddummy(EP2, false);
              if not(isempty(obj.Pind))
                EP1ind = [obj.Pind([i 1+rem(i,obj.NVert)],1); obj.Pind(i) + Ncube];
                EP2ind = [obj.Pind([i 1+rem(i,obj.NVert)],1); obj.Pind(1+rem(i,obj.NVert))];
                NewExtrudedFace1.Pind = EP1ind;
                NewExtrudedFace2.Pind = EP2ind;
              end
              ExtrudedFaces = [ExtrudedFaces; NewExtrudedFace1; NewExtrudedFace2]; %#ok 
           end
           EE = element3ddummy([obj.P; EP], ExtrudedFaces, false, [Pind; Pind + Ncube]);
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
        
    end
end

