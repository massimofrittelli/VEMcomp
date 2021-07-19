function [ KE,ME,volume,diameter] = matrices3Dlocal( E )

% Computes the local VEM stiffness and mass matrices of a 3D element.

% INPUTS:
% 1) an object E of class element3d

% OUTPUTS:
% 1) KE is the local stiffness matrix (size NVxNV)
% 2) ME is the local mass matrix (size NVxNV)
% 3) volume is the volume of the element
% 4) diameter is the diameter of the element

NV = E.Nvertices;

% %computing diameter and volume of the element
[volume, diameter] = barycenter3D(E); % TODO write function and convert to class method
% 
% %computing edges of the element
% edges = zeros(NV,2);
% for i=1:NV-1
%    edges(i,:) = PEP(i+1,:)-PEP(i,:); 
% end
% edges(NV,:) = PEP(1,:)-PEP(end,:);
% 
% %gradient of the non-constant scaled monomials 
% nabla1 = [1,0]/diameter;
% nabla2 = [0,1]/diameter;
% 
% %computing unit (outward ?) normals to each edge
% normals = -edges*[0,-1;1,0]';
% for i=1:NV
%    normals(i,:) = normals(i,:)/norm(normals(i,:)); 
% end
% 
% %computing outward normal derivatives of the non-constant scaled monomials
% %along the edges
% normders1 = normals * nabla1';
% normders2 = normals * nabla2';
% 
% %computing matrix B
% B = zeros(3,NV);
% B(1,:) = ones(1,NV)/NV;
% B(2,1) = (normders1(NV)*norm(edges(NV,:))+normders1(1)*norm(edges(1,:)))/2;
% B(3,1) = (normders2(NV)*norm(edges(NV,:))+normders2(1)*norm(edges(1,:)))/2;
% for i=2:NV
%    B(2,i) = (normders1(i-1)*norm(edges(i-1,:))+normders1(i)*norm(edges(i,:)))/2;
%    B(3,i) = (normders2(i-1)*norm(edges(i-1,:))+normders2(i)*norm(edges(i,:)))/2;
% end

%barycentric monomials
monomials = {@(x,y,z) 0*x + 1;
             @(x,y,z) (x-bary(1))/diameter;
             @(x,y,z) (y-bary(2))/diameter;
             @(x,y,z) (z-bary(3))/diameter};
         
%computing matrix D
D = zeros(NV,4);
D(:,1) = ones(NV,1);
for j=2:4
    for i=1:NV
        D(i,j) = monomials{j}(P(i,1),P(i,2),P(i,3));
    end
end

%COMPUTING LOCAL STIFFNESS MATRIX FROM B AND D (See Hitchhiker's)
G = B*D;
Gtilde = G;
Gtilde(1,:) = zeros(1,3);
Pnablastar = G\B;
Pnabla = D*Pnablastar;
KE = Pnablastar'*Gtilde*Pnablastar + diameter*(eye(NV)-Pnabla)'*(eye(NV)-Pnabla);

%computing matrix H         
H = zeros(4,4);
for i=1:4
    for j=1:4
        fun = @(x,y) monomials{i}(x,y).*monomials{j}(x,y);
        H(i,j) = quadrature_polyhedron(P, bary, fun);
    end
end

%COMPUTING LOCAL MASS MATRIX FROM H,B AND D (See Hitchhiker's)
C = H*Pnablastar;
P0 = Pnabla; %solo per k=1,2.
ME = C'*Pnablastar + volume*(eye(NV)-P0)'*(eye(NV)-P0);
end