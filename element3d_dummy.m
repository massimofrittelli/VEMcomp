classdef element3d_dummy
    %ELEMENT3D_DUMMY is a bare-bone class that models a convex 3D element:
    % 1) contains minimal information on the element
    % 2) allows for element translation/shifting
    % 3) allows for element plotting
    
    properties
        P(:,3) double % Vertexes
        Pind(:,1) double = [] % Indexes of vertexes
        Faces(:,1) element2d_dummy % Faces
        iscube(1,1) logical
    end
    
    properties (SetAccess = private)
        NVert(1,1) double
        NFaces(1,1) double
    end
    
    methods
        function obj = element3d_dummy(P,Faces,iscube,Pind)
            obj.P = P;
            obj.Faces = Faces;
            obj.iscube = iscube;
            if nargin >= 4
               obj.Pind = Pind; 
            end
            
            obj.NVert = length(P);
            obj.NFaces = length(Faces);
        end
        
        function P_ind_boundary = get_P_indexes_boundary(obj)
            P_ind_boundary = [];
            for i=1:length(obj.Faces)
                if obj.Faces(i).is_boundary
                   P_ind_boundary = [P_ind_boundary; obj.Faces(i).Pind]; %#ok
                end
            end
            P_ind_boundary = unique(P_ind_boundary);
        end
        
        function E = shiftElement(obj, v)
           newFaces = obj.Faces;
           for i=1:size(newFaces,1)
              newFaces(i) = shiftElement(obj.Faces(i), v); 
           end
           E =  element3d_dummy(obj.P+repmat(v,length(obj.P),1), newFaces, obj.iscube,obj.Pind);
        end
        
        function EE = extrude(obj, Ncube)
            % Extrudes the faces that need to be extruded
            if not(obj.iscube)
               error('Extrusion is only supported for cubic elements'); 
            end
            EE = [];
            I = eye(3);
            deltax = max(obj.P(:,1)) - min(obj.P(:,1));
            for i=1:6
                dir = obj.Faces(i).NormalDirection;
                fabsdir = abs(obj.Faces(i).P(1,dir));
                maxdir = max(obj.P(:,dir));
                mindir = min(obj.P(:,dir));
                maxabsdir = max(abs(obj.P(:,dir)));
                % discard faces that point toward the inside of the sphere
                % (3 out of 6)
                if maxdir*mindir > 0 && fabsdir < maxabsdir
                    continue
                end
                verse = sign(obj.Faces(i).P(1,dir) - mean(obj.P(:,dir)));
                % discard faces, based on normal direction, that are not to
                % be extruded
                if max(vecnorm((obj.Faces(i).P + I(dir,:)*deltax*verse)')) <=1
                   continue 
                end
                [~, indb] = ismember(obj.Faces(i).P,obj.P,'rows');
                obj.Faces(i).Pind = obj.Pind(indb);
                EE = [EE; extrude(obj.Faces(i), Ncube)]; %#ok
            end
        end
        
        function E = dummy2element(obj)
            if obj.iscube
               E = element3dcube(dummy2element(obj.Faces), obj.P, obj.Pind); 
            else
               E = element3d(dummy2element(obj.Faces), obj.P, mean(obj.P), obj.Pind); 
            end
        end
        
        function plot(obj, pert)
            if nargin == 1
               pert = 0; 
            end
           hold on
           for i=1:size(obj.Faces,1)
               PF = obj.Faces(i).P;
               % radii = vecnorm(PF')';
               fill3(PF(:,1), PF(:,2), PF(:,3), PF(:,1)*0 + 1 -pert - obj.iscube);
           end
        end
        
        function plotSolution(obj, sol)
           hold on
           solE = sol(obj.Pind);
           for i=1:size(obj.Faces,1)
               PF = obj.Faces(i).P;
               [tf,loc] = ismember(obj.P,PF,'rows');
               [~,p] = sort(loc(tf));
               indF = find(tf);
               indF = indF(p);
               if not(isempty(indF))
                    fill3(PF(:,1), PF(:,2), PF(:,3), solE(indF),'FaceColor', 'interp', 'EdgeColor', 'none');
               end
           end
        end
    end
end

