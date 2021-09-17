classdef element3ddummy
    %ELEMENT3DDUMMY is a bare-bone class that models a convex 3D element:
    % 1) contains minimal information on the element
    % 2) allows for element cutting in the unit sphere
    % 3) allows for element translation/shifting
    % 4) allows for element plotting
    
    properties
        P(:,3) % Vertexes
        Faces(:,1) element2ddummy % Faces
        iscube(1,1) logical
    end
    
    methods
        function obj = element3ddummy(P,Faces,iscube)
            obj.P = P;
            obj.Faces = Faces;
            obj.iscube = iscube;
        end
        
        function E = shiftElement(obj, v)
           newFaces = obj.Faces;
           for i=1:size(newFaces,1)
              newFaces(i) = shiftElement(obj.Faces(i), v); 
           end
           E =  element3ddummy(obj.P+repmat(v,length(obj.P),1), newFaces, obj.iscube);
        end
        
        function E = dummy2element(obj)
            if obj.iscube
               E = element3dcube(dummy2element(obj.Faces), obj.P); 
            else
               E = element3d(dummy2element(obj.Faces), obj.P, mean(obj.P)); 
            end
        end
        
        function obj = cut(obj)
            % Works only for cutting convex elements to fit in the sphere
            
            newFaces = [];
            for i=1:length(obj.Faces)
                newFaces = [newFaces; cut(obj.Faces(i))]; %#ok
            end
            obj.Faces = newFaces;
            
            for i=1:size(obj.P,1)
               if norm(obj.P(i,:)) > 1
                   obj.iscube = false;
                   obj.P(i,:) = obj.P(i,:)/norm(obj.P(i,:));
               end
            end
        end
        
        function plot(obj)
           for i=1:size(obj.Faces,1)
               PF = obj.Faces(i).P;
               radii = vecnorm(PF')';
               fill3(PF(:,1), PF(:,2), PF(:,3), radii.^3);
           end
        end
    end
end

