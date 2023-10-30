% Generating mesh
fun = @(P) P(:,1).^2 + P(:,2).^2 + P(:,3).^2 -1;
range = [-1,1; -1,1; -1,1];
tol = 1e-3;
Nx = 30;
[P, h, BulkElements, SurfaceElements, ElementsPlot] = ...
    generate_mesh_3d(fun, range, Nx, tol, -0.3);

% Assembling matrices
[K, M, C, KS, MS, CS, R] = assembly_3d(P, BulkElements, SurfaceElements);

% Solving PDE
% Neumann
% f = {@(P) 4*(3-5*(P(:,1).^2 + P(:,2).^2 + P(:,3).^2)) ...
%    + (1- (P(:,1).^2 + P(:,2).^2 + P(:,3).^2)).^2};
% Dirichlet
f = {@(P) 7 - (P(:,1).^2 + P(:,2).^2 + P(:,3).^2)};
D = 1;
alpha = 1;
u = solver_elliptic_bulk(1, D, alpha, f, P, M, K, R, 'dir');

% Computing relative L2 error
% Neumann
% esol = @(P) (1- (P(:,1).^2 + P(:,2).^2 + P(:,3).^2)).^2;
% Dirichlet
esol = @(P) 1- (P(:,1).^2 + P(:,2).^2 + P(:,3).^2);
u_exact = esol(P);
pointwise_error = u_exact - u;
normsol = sqrt(u_exact'*M*u_exact);
normerror = sqrt(pointwise_error'*M*pointwise_error);
relative_err = normerror/normsol;

% Plotting numerical solution
figure
set(gcf,'color','white')
for i=1:length(ElementsPlot)
   plot(ElementsPlot(i), u(ElementsPlot(i).Pind)); 
   hold on
end
colormap jet
axis equal tight
xlabel('x')
ylabel('y','rot',0)
zlabel('z','rot',0)
title('u')
set(gca,'FontSize',18)
colorbar

figure
set(gcf,'color','white')
for i=1:length(ElementsPlot)
   plot(ElementsPlot(i), esol(P(ElementsPlot(i).Pind, :))); 
   hold on
end
colormap jet
axis equal tight
xlabel('x')
ylabel('y','rot',0)
zlabel('z','rot',0)
set(gca,'FontSize',18)
title('u_{exact}')
colorbar