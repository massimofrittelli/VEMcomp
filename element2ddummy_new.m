classdef element2ddummy_new
    %ELEMENT2DDUMMY is a bare-bone class that models a convex 2D element:
    % 1) contains minimal information on the element
    % 2) allows for element extrusion in the unit sphere
    % 3) allows for element translation/shifting
    
    properties
        Pind(:,1) double = [] % Indexes of vertexes
        boundary(1,1) logical
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
        function obj = element2ddummy_new(P,issquare,boundary,Pind)     
            obj.P = P;
            obj.NVert = size(P,1);
            obj.issquare = issquare;
            if nargin >= 3
                obj.boundary = boundary;
            end
            if nargin >= 4
                obj.Pind = Pind;
            end
        end
        
        function ND = get.NormalDirection(obj)
           if not(obj.issquare)
              error('Normal direction is supported square elements only') 
           end
           direction = cross(obj.P(3,:)-obj.P(2,:) , obj.P(2,:)-obj.P(1,:));
           [~, ND] = max(abs(direction));
        end
        
        function E = shiftElement(obj, v)
           E =  element2ddummy_new(obj.P+repmat(v,length(obj.P),1),obj.issquare,false,obj.Pind);
        end
        
        
        function EE = extrude(obj, Ncube)
           EP = zeros(size(obj.P));
           actuallyExtruded = logical(size(obj.P,1));
           for i=1:obj.NVert
               [EP(i,:), actuallyExtruded(i)] = extrude_node(obj.P(i,:));
           end
           extruded_ind = obj.Pind + Ncube;
           extruded_ind(not(actuallyExtruded)) = obj.Pind(not(actuallyExtruded));
           ExtrudedFaces = obj;
           % CREATE TRIANGULAR FACES LYING ON THE SURFACE
           if norm(EP(1,:) + EP(3,:)) >= norm(EP(2,:) + EP(4,:))
               NewExtrudedFace1 = element2ddummy_new(EP([1 2 3],:), false, true);
               NewExtrudedFace1.Pind = extruded_ind([1 2 3],1);
               NewExtrudedFace2 = element2ddummy_new(EP([1 3 4],:), false, true);
               NewExtrudedFace2.Pind = extruded_ind([1 3 4],1);
           else
               NewExtrudedFace1 = element2ddummy_new(EP([1 2 4],:), false, true);
               NewExtrudedFace1.Pind = extruded_ind([1 2 4],1);
               NewExtrudedFace2 = element2ddummy_new(EP([2 3 4],:), false, true);
               NewExtrudedFace2.Pind = extruded_ind([2 3 4],1);
           end
               ExtrudedFaces = [ExtrudedFaces; NewExtrudedFace1; NewExtrudedFace2];
           % CREATE LATERAL SQUARE FACES
           for i=1:obj.NVert
              LP = unique([obj.P([i 1+rem(i,obj.NVert)],:); EP([1+rem(i,obj.NVert) i],:)],'rows','stable');
              NewExtrudedFace = element2ddummy_new(LP, false, false);
              EPind = unique([obj.Pind([i 1+rem(i,obj.NVert)],1); extruded_ind([1+rem(i,obj.NVert) i],1)],'stable');
              NewExtrudedFace.Pind = EPind;
              ExtrudedFaces = [ExtrudedFaces; NewExtrudedFace]; %#ok 
           end
           EE = element3ddummy_new([obj.P; EP(actuallyExtruded,:)], ExtrudedFaces, false, [obj.Pind; extruded_ind(actuallyExtruded)], extruded_ind);
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

function [P, wasActuallyExtruded] = extrude_node(P)
    tol = 1e-8;
    % TYPE 1: LESS UNIFORM, BUT INTERNAL FACES ARE GUARANTEED TO BE FLAT
    n = signtol(P);
    alpha = (- P*n' + sqrt((P*n')^2 -norm(n)^2*(norm(P)^2-1)))/norm(n)^2;
    wasActuallyExtruded = (abs(alpha) > tol);
    P = P + wasActuallyExtruded*alpha*n;
    % TYPE 2: MORE UNIFORM, BUT INTERNAL FACES ARE NOT GUARANTEED TO BE
    % FLAT
%     wasActuallyExtruded = (1 - norm(P) > tol);
%     P = P/norm(P);
end

function x = signtol(x)
    tol = 1e-10;
    x = sign(x).*(abs(x) > tol);
end
