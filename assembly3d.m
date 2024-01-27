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
% - K,M,C: stiffness, mass, and consistency matrices in the bulk
% - KS,MS,CMS: stiffness, mass, and consistency matrices on the surface
% - R: reduction matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [K, C, M, KS, CS, MS, R] = assembly3d(P, BulkElements, SurfElements)

Nbulk = length(P);
boundarynodes = unique(SurfElements(:));
Nsurf = length(boundarynodes);

if isempty(BulkElements) % Triangulated surface mesh
    K = []; M = []; C = []; R = speye(Nsurf);
    [KS,MS] = matrices(P,SurfElements);
    CS = MS;
    return
end

K = spalloc(Nbulk,Nbulk,27*Nbulk); % Stiffness matrix in the bulk
M = spalloc(Nbulk,Nbulk,27*Nbulk); % Mass matrix in the bulk
C = spalloc(Nbulk,Nbulk,27*Nbulk); % Mass matrix in the bulk

KS = spalloc(Nbulk,Nbulk,9*Nsurf); % Stiffness matrix on the surf
MS = spalloc(Nbulk,Nbulk,9*Nsurf); % Mass matrix on the surf
CS = spalloc(Nbulk,Nbulk,9*Nsurf); % Consistency matrix on the surf

% find first cubic element in mesh (they are all equal)
for i=1:length(BulkElements)
    if BulkElements(i).is_cube
        Cube = copyElement3d(BulkElements(i));
        MC = getLocalMatrices(Cube).M;
        KC = getLocalMatrices(Cube).K;
        CC = getLocalMatrices(Cube).C;
        % An element3dcube is not supposed to have boundary faces
        break
    end
end

% MATRIX ASSEMBLY
for i=1:length(BulkElements) % For each bulk element
    ElementDummy = BulkElements(i);
    eind = ElementDummy.Pind;
    if ElementDummy.is_cube
        M(eind, eind) = M(eind, eind) + MC; %#ok
        C(eind, eind) = C(eind, eind) + CC; %#ok
        K(eind, eind) = K(eind, eind) + KC; %#ok
    else
        Element = copyElement3d(ElementDummy);
        try
        M(eind, eind) = M(eind, eind) + getLocalMatrices(Element).M; %#ok
        catch
        print('ops');
        end
        C(eind, eind) = C(eind, eind) + getLocalMatrices(Element).C; %#ok
        K(eind, eind) = K(eind, eind) + getLocalMatrices(Element).K; %#ok
        for j=1:Element.NFaces
            Face = copyElement2d(Element.Faces(j));
            if Face.is_boundary
                eind_boundary = Face.Pind;
                MS(eind_boundary, eind_boundary) = MS(eind_boundary, eind_boundary) + getLocalMatrices(Face).M; %#ok
                KS(eind_boundary, eind_boundary) = KS(eind_boundary, eind_boundary) + getLocalMatrices(Face).K; %#ok
                CS(eind_boundary, eind_boundary) = CS(eind_boundary, eind_boundary) + getLocalMatrices(Face).C; %#ok
            end
        end
    end
end

MS = MS(boundarynodes,boundarynodes);
CS = CS(boundarynodes,boundarynodes);
KS = KS(boundarynodes,boundarynodes);
            
R = spalloc(Nbulk, Nsurf, Nsurf);
R(boundarynodes,:) = speye(Nsurf);

        
end