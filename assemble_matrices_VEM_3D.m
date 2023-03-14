%
% DESCRIPTION - Assembles bulk and surface VEM meshes given a polyhedral
% mesh
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INPUTS:
%
%  - P: array of nodes
%  - Elements: all elements of the bulk
%  - EGamma: all elements of the surface
%
% OUTPUTS:
%
% - K,M: stiffness and mass matrices in the bulk
% - KS,MS,CMS: stiffness, mass, and consistency matrices on the surface
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [K, M, KS, MS, CMS] = assemble_matrices_VEM_3D(P, EGamma, Elements)

Nbulk = length(P);
Nsurf = max(max(EGamma));

K = spalloc(Nbulk,Nbulk,57*Nbulk); % Stiffness matrix in the bulk
M = spalloc(Nbulk,Nbulk,57*Nbulk); % Mass matrix in the bulk

KS = spalloc(Nsurf,Nsurf,9*Nsurf); % Stiffness matrix on the surf
MS = spalloc(Nsurf,Nsurf,9*Nsurf); % Mass matrix on the surf
CMS = spalloc(Nsurf,Nsurf,9*Nsurf); % Consistency matrix on the surf

% find first cubic element in mesh (they are all equal)
for i=1:length(Elements)
    if Elements(i).iscube
        Cube = dummy2element(Elements(i));
        MC = Cube.M;
        KC = Cube.K;
        % An element3dcube is not supposed to have boundary faces
        break
    end
end

% MATRIX ASSEMBLY
for i=1:length(Elements) % For each bulk element
    ElementDummy = Elements(i);
    eind = ElementDummy.Pind;
    if ElementDummy.iscube
        M(eind, eind) = M(eind, eind) + MC; %#ok
        K(eind, eind) = K(eind, eind) + KC; %#ok
    else
        Element = dummy2element(ElementDummy);
        M(eind, eind) = M(eind, eind) + Element.M; %#ok
        K(eind, eind) = K(eind, eind) + Element.K; %#ok
        for j=1:Element.NFaces
            Face = Element.Faces(j);
            if Face.is_boundary
                eind_boundary = Face.Pind;
                MS(eind_boundary, eind_boundary) = MS(eind_boundary, eind_boundary) + Face.M; %#ok
                KS(eind_boundary, eind_boundary) = KS(eind_boundary, eind_boundary) + Face.K; %#ok
                CMS(eind_boundary, eind_boundary) = CMS(eind_boundary, eind_boundary) + Face.CM; %#ok
            end
        end
    end
end
            
        
end