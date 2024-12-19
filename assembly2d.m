%
% DESCRIPTION - Assembles bulk and surface VEM meshes given a polygonal
% mesh
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INPUTS:
%
%  - P: array of nodes
%  - BulkElements: all elements of the bulk
%  - SurfaceElements: all elements of the surface
%
% OUTPUTS:
%
% - K,C,M: stiffness, consistency and mass matrices in the bulk
% - KS,MS: stiffness and mass matrices on the surface
% - R: reduction matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [K,C,M, KS, MS, R] = assembly2d(P, BulkElements, SurfElements)

Nbulk = length(P);
boundarynodes = unique(SurfElements(:));
Nsurf = length(boundarynodes);

% find first square element in mesh (they are all equal)
for i=1:length(BulkElements)
    if BulkElements(i).is_square
        Square = getLocalMatrices(copyElement2d(BulkElements(i)));
        MSq = Square.M(:);
        CSq = Square.C(:);
        KSq = Square.K(:);
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
Nsurf_repeated = 4*Nsurf;
sii = zeros(Nsurf_repeated,1);
sjj = zeros(Nsurf_repeated,1);
svk = zeros(Nsurf_repeated,1);
svm = zeros(Nsurf_repeated,1);

% MATRIX ASSEMBLY IN THE BULK
bulk_counter = 0;
for i=1:length(BulkElements) % For each bulk element
    Element = copyElement2d(BulkElements(i));
    eind = Element.Pind;
    oind = ones(Element.NVert,1);
    ii(bulk_counter+1: bulk_counter+Element.NVert^2) = kron(oind,eind);
    jj(bulk_counter+1: bulk_counter+Element.NVert^2) = kron(eind,oind);
    if Element.is_square
        vm(bulk_counter+1: bulk_counter+Element.NVert^2) = MSq;
        vc(bulk_counter+1: bulk_counter+Element.NVert^2) = CSq;
        vk(bulk_counter+1: bulk_counter+Element.NVert^2) = KSq;
    else
        Element = getLocalMatrices(Element);
        vm(bulk_counter+1: bulk_counter+Element.NVert^2) = Element.M(:);
        vc(bulk_counter+1: bulk_counter+Element.NVert^2) = Element.C(:);
        vk(bulk_counter+1: bulk_counter+Element.NVert^2) = Element.K(:);
    end
    bulk_counter = bulk_counter + Element.NVert^2;
end

% MATRIX ASSEMBLY ON THE SURFACE
surf_counter = 0;
for i=1:length(SurfElements)
    eind = SurfElements(i,:)';
    oind = [1;1];
    sii(surf_counter+1: surf_counter+4) = kron(eind,oind);
    sjj(surf_counter+1: surf_counter+4) = kron(oind,eind);
    element_length = norm(P(eind(1),:) - P(eind(2),:));
    svm(surf_counter+1: surf_counter+4) = [2 1 1 2]'/6*element_length;
    svk(surf_counter+1: surf_counter+4) = [1 -1 -1 1]'/element_length;
    surf_counter = surf_counter + 4;
end

MS = sparse(sii,sjj,svm);
KS = sparse(sii,sjj,svk);
MS = MS(boundarynodes,boundarynodes);
KS = KS(boundarynodes,boundarynodes);
            
R = spalloc(Nbulk, Nsurf, Nsurf);
R(boundarynodes,:) = speye(Nsurf);

M = sparse(ii,jj, vm);
C = sparse(ii,jj, vc);
K = sparse(ii,jj, vk);

end