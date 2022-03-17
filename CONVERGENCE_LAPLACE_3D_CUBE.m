% DESCRIPTION - Tests the convergence for the Laplace equation on the 3D
% cube with homogeneous Neumann boundary conditions with VEM on a cubic
% mesh.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NX = [6,11,16,21,26,31,36,41,46];
Nit = length(NX);

L2errors = zeros(1,Nit);
L2errors_rel = zeros(1,Nit);
meshsizes = zeros(1,Nit);
L2rates = zeros(1,Nit-1);
L2rates_rel = zeros(1,Nit-1);

for i=1:Nit-1
   [L2errors(i), L2errors_rel(i), meshsizes(i)] = laplace3dcube(NX(i), false); 
end
[L2errors(Nit), L2errors_rel(Nit), meshsizes(Nit)] = laplace3dcube(NX(Nit), true);

for i=1:Nit-1
   L2rates(i) = log(L2errors(i+1)/L2errors(i)) / log(meshsizes(i+1)/meshsizes(i));
   L2rates_rel(i) = log(L2errors_rel(i+1)/L2errors_rel(i)) / log(meshsizes(i+1)/meshsizes(i));
end