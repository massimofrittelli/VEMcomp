function [K,M] = assemble_matrices_VEM_2D(P, Elements)
%ASSEMBLE_MATRICES_VEM_2D Summary of this function goes here
%   Detailed explanation goes here

Nbulk = length(P);

K = spalloc(Nbulk,Nbulk,57*Nbulk); % Stiffness matrix in the 2D bulk
M = spalloc(Nbulk,Nbulk,57*Nbulk); % Mass matrix in the 2D bulk

% find first square element in mesh (they are all equal)
for i=1:length(Elements)
    if Elements(i).is_square
        Square = dummy2element(Elements(i));
        MS = Square.M;
        KS = Square.K;
        % An element2dsquare is not supposed to have boundary faces
        break
    end
end

% MATRIX ASSEMBLY
for i=1:length(Elements) % For each bulk element
    ElementDummy = Elements(i);
    eind = ElementDummy.Pind;
    if ElementDummy.is_square
        M(eind, eind) = M(eind, eind) + MS; %#ok
        K(eind, eind) = K(eind, eind) + KS; %#ok
    else
        Element = dummy2element(ElementDummy);
        M(eind, eind) = M(eind, eind) + Element.M; %#ok
        K(eind, eind) = K(eind, eind) + Element.K; %#ok
    end
end

end