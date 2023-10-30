% Generating mesh
fun = @(P) P(:,1).^2 + P(:,2).^2 -1;
range = [-1,1; -1,1];
tol = 1e-3;
Nx = 40;
[P, h, BulkElements, SurfaceElements] = generate_mesh_2d(fun, range, Nx, tol);

% Assembling matrices
[K,C,M,KS,MS,R] = assembly_2d(P, BulkElements, SurfaceElements);

% Solving PDE
% Neumann
f = {@(P) 8*(1-2*(P(:,1).^2 + P(:,2).^2)) + (1- (P(:,1).^2 + P(:,2).^2)).^2};
% Dirichlet 
% f = {@(P) 5 - P(:,1).^2 - P(:,2).^2};
D = 1;
alpha = 1;
u = solver_elliptic_bulk(1, D, alpha, f, P, M, K, R, 'neu');

% Computing relative L2 error
% Neumann
esol = @(P) (1- (P(:,1).^2 + P(:,2).^2)).^2;
% Dirichlet
% esol = @(P) 1- (P(:,1).^2 + P(:,2).^2);
u_exact = esol(P);
pointwise_error = u_exact - u;
normsol = sqrt(u_exact'*M*u_exact);
normerror = sqrt(pointwise_error'*M*pointwise_error);
relative_err = normerror/normsol;

% Plotting numerical solution
figure
set(gcf,'color','white')
for i=1:length(BulkElements)
   plot(BulkElements(i), u(BulkElements(i).Pind)); 
   hold on
end
colormap jet
view(2)
axis equal tight
xlabel('x')
ylabel('y','rot',0)
title('u_h')
set(gca,'FontSize',18)
colorbar

figure
set(gcf,'color','white')
for i=1:length(BulkElements)
   plot(BulkElements(i), esol(P(BulkElements(i).Pind, :))); 
   hold on
end
colormap jet
view(2)
axis equal tight
xlabel('x')
ylabel('y','rot',0)
set(gca,'FontSize',18)
title('u_{exact}')
colorbar