classdef element3ddummy
    %ELEMENT3DDUMMY is a bare-bone class that models a convex 3D element:
    % 1) contains minimal information on the element
    % 2) allows for element cutting in the unit sphere
    % 3) allows for element plotting
    
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
           figure
           hold on
           for i=1:size(obj.Faces,1)
               PF = obj.Faces(i).P;
               fill3(PF(:,1), PF(:,2), PF(:,3), 0*PF(:,1));
           end
           view(3)
           axis equal
           xlabel('x')
           ylabel('y')
           zlabel('z','rot',0)
        end
    end
end

