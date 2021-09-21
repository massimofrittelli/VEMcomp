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
        
        function EE = extrude(obj, Pind, Ncube)
            % Extrudes the faces that contain the point OP as vertex in the
            % directions specified by the vector 'directions'
            if not(obj.iscube)
               error('Extrusion is only supported for cubic elements'); 
            end
            obj.Pind = Pind; % Indexes of vertexes of obj
            EE = [];
            I = eye(3);
            deltax = max(obj.P(:,1)) - min(obj.P(:,1));
            for i=1:6
                dir = obj.Faces(i).NormalDirection;
                fabsdir = abs(obj.Faces(i).P(1,dir));
                maxdir = max(obj.P(:,dir));
                mindir = min(obj.P(:,dir));
                maxabsdir = max(abs(obj.P(:,dir)));
                if maxdir*mindir > 0 && fabsdir < maxabsdir
                    continue
                end
                verse = sign(obj.Faces(i).P(1,dir) - mean(obj.P(:,dir)));
                if max(vecnorm((obj.Faces(i).P + I(dir,:)*deltax*verse)')) <=1
                   continue 
                end
                [~, indb] = ismember(obj.Faces(i).P,obj.P,'rows');
                EE = [EE; extrude(obj.Faces(i), Pind(indb), Ncube)]; %#ok
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
               % radii = vecnorm(PF')';
               fill3(PF(:,1), PF(:,2), PF(:,3), PF(:,1)*0 + 1 - obj.iscube);
           end
        end
    end
end

