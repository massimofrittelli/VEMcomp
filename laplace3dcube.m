% DESCRIPTION - Solves Laplace equation on the unit cube with homogeneous
% Neumann boundary conditions with VEM on a cubic mesh.
%
%   - Delta(u) + u = (3*pi^2+1)*cos(pi*x)*cos(pi*y)*cos(pi*z)
%     Nabla(u) \cdot n = 0
%
% Exact solution: u(x,y,z) = cos(pi*x)*cos(pi*y)*cos(pi*z)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INPUTS:
%
% - Nx: number of nodes along each direction
% - plotflag: indicates whether generate plots
%
% OUTPUTS:
%
% - l2err: absolute error in L2 norm
% - l2err_rel: relative error in L2 norm
% - h: meshsize
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [l2err,l2err_rel, h] = laplace3dcube(Nx,plotflag)

[P, h, K, M, ~, ~, ~, ~] = generate_mesh_cube(Nx);

f = @(x,y,z) (3*pi^2+1)*cos(pi*x).*cos(pi*y).*cos(pi*z);
esol = @(x,y,z) cos(pi*x).*cos(pi*y).*cos(pi*z);

% Solving linear system (K+M)u = Mf
tic
u = (K+M)\(M*f(P(:,1),P(:,2),P(:,3)));
es = esol(P(:,1),P(:,2),P(:,3));
toc

% Computing L2 error
err = es - u;
l2err = sqrt(err'*M*err);
l2solnorm = sqrt(es'*M*es);
l2err_rel = l2err/l2solnorm;

if plotflag
    
% Plotting Exact solution
figure
set(gcf,'Color','white')
scatter3(P(:,1), P(:,2), P(:,3), 30, esol(P(:,1), P(:,2), P(:,3)),'filled')
set(gca,'FontSize',18)
xlabel('x')
ylabel('y')
zlabel('z','rot',0)
title('Exact solution')
colorbar

% Plotting Numerical solution
figure
set(gcf,'Color','white')
scatter3(P(:,1), P(:,2), P(:,3), 30, u,'filled')
set(gca,'FontSize',18)
xlabel('x')
ylabel('y')
zlabel('z','rot',0)
title('Numerical solution')
colorbar

end

end