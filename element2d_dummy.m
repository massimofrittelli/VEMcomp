classdef element2d_dummy
    %ELEMENT2D_DUMMY is a bare-bone class that models a convex 2D element:
    % 1) contains minimal information on the element
    % 2) allows for element extrusion in the unit sphere
    % 3) allows for element translation/shifting
    
    properties
        P(:,3) double % Vertexes
        Pind(:,1) double = [] % Indexes of vertexes
        to_plot logical
    end
    
    properties (SetAccess = private)
        NVert(1,1) double % Number of vertexes
        is_square(1,1) logical
        is_boundary(1,1) logical
    end
    
    properties (Dependent)
       NormalDirection(3,1) double 
    end
    
    methods
        function obj = element2d_dummy(P,is_square,is_boundary,Pind)     
            obj.P = P;
            obj.NVert = size(P,1);
            obj.is_square = is_square;
            if nargin >= 3
                obj.is_boundary = is_boundary;
            end
            if nargin >= 4
                obj.Pind = Pind;
            end
        end
        
        function ND = get.NormalDirection(obj)
           if not(obj.is_square)
              error('Normal direction is supported square elements only') 
           end
           direction = cross(obj.P(3,:)-obj.P(2,:) , obj.P(2,:)-obj.P(1,:));
           [~, ND] = max(abs(direction));
        end
        
        function E = shiftElement(obj, v)
           E =  element2d_dummy(obj.P+repmat(v,length(obj.P),1),obj.is_square,false,obj.Pind);
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
           
           % Re-ordering vertexes of square face of extruded element, in
           % such a wat that the normal is outward
           ExtrudedFaces.P = flipud(ExtrudedFaces.P);
           ExtrudedFaces.Pind = flipud(ExtrudedFaces.Pind);
           
           % CREATE TRIANGULAR FACES LYING ON THE SURFACE
           if norm(EP(1,:) + EP(3,:)) >= norm(EP(2,:) + EP(4,:))
               NewExtrudedFace1 = element2d_dummy(EP([1 2 3],:), false, true);
               NewExtrudedFace1.Pind = extruded_ind([1 2 3],1);
               NewExtrudedFace2 = element2d_dummy(EP([1 3 4],:), false, true);
               NewExtrudedFace2.Pind = extruded_ind([1 3 4],1);
           else
               NewExtrudedFace1 = element2d_dummy(EP([1 2 4],:), false, true);
               NewExtrudedFace1.Pind = extruded_ind([1 2 4],1);
               NewExtrudedFace2 = element2d_dummy(EP([2 3 4],:), false, true);
               NewExtrudedFace2.Pind = extruded_ind([2 3 4],1);
           end
               ExtrudedFaces = [ExtrudedFaces; NewExtrudedFace1; NewExtrudedFace2];
           % CREATE LATERAL SQUARE FACES
           for i=1:obj.NVert
              LP = unique([obj.P([i 1+rem(i,obj.NVert)],:); EP([1+rem(i,obj.NVert) i],:)],'rows','stable');
              NewExtrudedFace = element2d_dummy(LP, false, false);
              EPind = unique([obj.Pind([i 1+rem(i,obj.NVert)],1); extruded_ind([1+rem(i,obj.NVert) i],1)],'stable');
              NewExtrudedFace.Pind = EPind;
              ExtrudedFaces = [ExtrudedFaces; NewExtrudedFace]; %#ok 
           end
           EE = element3d_dummy([obj.P; EP(actuallyExtruded,:)], ExtrudedFaces, false, [obj.Pind; extruded_ind(actuallyExtruded)]);
        end
        
        function Earray = dummy2element(EDarray)
            for ei = 1:length(EDarray)
                if EDarray(ei).is_square
                    Earray(ei) = element2dsquare(EDarray(ei).P, EDarray(ei).Pind, EDarray(ei).is_boundary); %#ok
                else
                    Earray(ei) = element2d(EDarray(ei).P, mean(EDarray(ei).P), EDarray(ei).Pind, EDarray(ei).is_boundary); %#ok
                end
            end
        end
        
        function plot(obj, faceColor, edgeColor)
            switch nargin
                case 1
                    fill3(obj.P(:,1), obj.P(:,2), obj.P(:,3), obj.P(:,1)*0 + 1 - obj.is_square);
                case 2
                    fill3(obj.P(:,1), obj.P(:,2), obj.P(:,3), faceColor);
                case 3
                    fill3(obj.P(:,1), obj.P(:,2), obj.P(:,3), faceColor, 'EdgeColor', edgeColor);
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
