classdef element3ddummy
    %ELEMENT3DDUMMY is a bare-bone class that models a convex 3D element:
    % 1) contains minimal information on the element
    % 2) allows for element translation/shifting
    % 3) allows for element plotting
    
    properties
        P(:,3) double % Vertexes
        Pind(:,1) double = [] % Indexes of vertexes
        Faces(:,1) element2ddummy % Faces
        iscube(1,1) logical
    end
    
    methods
        function obj = element3ddummy(P,Faces,iscube,Pind)
            obj.P = P;
            obj.Faces = Faces;
            obj.iscube = iscube;
            if nargin == 4
               obj.Pind = Pind; 
            end
        end
        
        function E = shiftElement(obj, v)
           newFaces = obj.Faces;
           for i=1:size(newFaces,1)
              newFaces(i) = shiftElement(obj.Faces(i), v); 
           end
           E =  element3ddummy(obj.P+repmat(v,length(obj.P),1), newFaces, obj.iscube);
        end
        
        function EE = extrude(obj, Pind, directions, verses, Ncube)
            % Extrudes the faces that contain the point OP as vertex in the
            % directions specified by the vector 'directions'
            if not(obj.iscube)
               error('Extrusion is only supported for cubic elements'); 
            end
            obj.Pind = Pind; % Indexes of vertexes of obj
            EE = [];
            for i=1:6
               dirindex = find(ismember(directions, obj.Faces(i).NormalDirection));
               if isempty(dirindex)
                  continue 
               end
               dir = directions(dirindex);
               verse = verses(dirindex);
               if verse*(obj.Faces(i).P(1,dir) - mean(obj.P(:,dir))) > 0
                   [~, indb] = ismember(obj.Faces(i).P,obj.P,'rows');
                   EE = [EE; extrude(obj.Faces(i), Pind(indb), Ncube)]; %#ok
               end
            end
        end
        
        function E = dummy2element(obj)
            if obj.iscube
               E = element3dcube(dummy2element(obj.Faces), obj.P); 
            else
               E = element3d(dummy2element(obj.Faces), obj.P, mean(obj.P)); 
            end
        end
        
        function plot(obj)
           hold on
           for i=1:size(obj.Faces,1)
               PF = obj.Faces(i).P;
               radii = vecnorm(PF')';
               fill3(PF(:,1), PF(:,2), PF(:,3), radii.^3);
           end
        end
    end
end

