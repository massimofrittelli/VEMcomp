close all
clearvars

% Generating mesh
fun = @(P) P(:,1).^2 + P(:,2).^2 -1;
range = [-1,1; -1,1];
tol = 1e-6;
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
L2_relative_err = compute_error(C,MS,u,u_exact,[],[]);

% Plotting numerical solution
figure
set(gcf,'color','white')
plot_bulk_2d(BulkElements, u, '$u$')

% PLOTTING EXACT SOLUTION
figure
set(gcf,'color','white')
plot_bulk_2d(BulkElements, esol(P), '$u_{exact}$')