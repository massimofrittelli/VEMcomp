% % Generating mesh
fun = @(P) P(:,1).^2 + P(:,2).^2 + P(:,3).^2 -1;
range = [-1,1; -1,1; -1,1]*1.04;
tol = 1e-10;
Nx = 10;
[P, h, BulkElements, SurfElements, ElementsPlot] = ...
     generate_mesh3d(fun, range, Nx, tol, -0.3);
% 
% % Assembling matrices
[K, M, C, KS, MS, CS, R] = assembly3d(P, BulkElements, SurfElements);

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
[L2_relative_err] = compute_error(C,[],u,u_exact,[],[]);

% Plotting numerical solution
figure, set(gcf,'color','white')
plot_bulk3d(ElementsPlot, u, '$u$');

% Plotting exact solution
figure, set(gcf,'color','white')
plot_bulk3d(ElementsPlot, esol(P), '$u_{exact}$');