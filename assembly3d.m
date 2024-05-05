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

% find first cubic element in mesh (they are all equal)
for i=1:length(BulkElements)
    if BulkElements(i).is_cube
        Cube = getLocalMatrices(copyElement3d(BulkElements(i)));
        MC = Cube.M(:);
        KC = Cube.K(:);
        CC = Cube.C(:);
        % An element3dcube is not supposed to have boundary faces
        break
    end
end

% allocate vectors for sparse matrix representation
Nbulk_repeated = 0;
for i=1:length(BulkElements)
    Nbulk_repeated = Nbulk_repeated + BulkElements(i).NVert;
end
ii = zeros(Nbulk_repeated,1);
jj = zeros(Nbulk_repeated,1);
vk = zeros(Nbulk_repeated,1);
vm = zeros(Nbulk_repeated,1);
vc = zeros(Nbulk_repeated,1);
Nsurf_repeated = size(SurfElements,1)*3;
sii = zeros(Nsurf_repeated,1);
sjj = zeros(Nsurf_repeated,1);
svk = zeros(Nsurf_repeated,1);
svm = zeros(Nsurf_repeated,1);
svc = zeros(Nsurf_repeated,1); 


% MATRIX ASSEMBLY
bulk_counter = 0;
surf_counter = 0;
for i=1:length(BulkElements) % For each bulk element
    E = BulkElements(i);
    eind = E.Pind;
    oind = ones(E.NVert,1);
    ii(bulk_counter+1: bulk_counter+E.NVert^2) = kron(oind,eind);
    jj(bulk_counter+1: bulk_counter+E.NVert^2) = kron(eind,oind);
    if E.is_cube
        vm(bulk_counter+1: bulk_counter+E.NVert^2) = MC;
        vc(bulk_counter+1: bulk_counter+E.NVert^2) = CC;
        vk(bulk_counter+1: bulk_counter+E.NVert^2) = KC;
    else
        Element = getLocalMatrices(copyElement3d(E));
        vm(bulk_counter+1: bulk_counter+E.NVert^2) = Element.M(:);
        vc(bulk_counter+1: bulk_counter+E.NVert^2) = Element.C(:);
        vk(bulk_counter+1: bulk_counter+E.NVert^2) = Element.K(:);
        for j=1:Element.NFaces
            Face = copyElement2d(Element.Faces(j));
            if Face.is_boundary
                eind_boundary = Face.Pind;
                oind_boundary = ones(Face.NVert,1);
                sii(surf_counter+1: surf_counter+Face.NVert^2) = kron(oind_boundary,eind_boundary);
                sjj(surf_counter+1: surf_counter+Face.NVert^2) = kron(eind_boundary,oind_boundary);
                Face = getLocalMatrices(Face);
                svm(surf_counter+1: surf_counter+Face.NVert^2) = Face.M(:);
                svc(surf_counter+1: surf_counter+Face.NVert^2) = Face.C(:);
                svk(surf_counter+1: surf_counter+Face.NVert^2) = Face.K(:);
                surf_counter = surf_counter + Face.NVert^2;
            end
        end
    end
    bulk_counter = bulk_counter + E.NVert^2;
end

MS = sparse(sii,sjj, svm);
CS = sparse(sii,sjj, svc);
KS = sparse(sii,sjj, svk);
MS = MS(boundarynodes,boundarynodes);
CS = CS(boundarynodes,boundarynodes);
KS = KS(boundarynodes,boundarynodes);
            
R = spalloc(Nbulk, Nsurf, Nsurf);
R(boundarynodes,:) = speye(Nsurf);

M = sparse(ii,jj, vm);
C = sparse(ii,jj, vc);
K = sparse(ii,jj, vk);
        
end