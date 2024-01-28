close all
clearvars

% STEP 1: Generate mesh
level_fun = @(P) P(:,1).^2 + P(:,2).^2 -1;
range = [-1,1; -1,1]; Nx = 40; tol = 1e-6;
[P, h, BulkElements, SurfElements] = generate_mesh2d(level_fun, range, Nx, tol);

% Step 2: Matrix assembly
[K,C,M,KS,MS,R] = assembly2d(P, BulkElements, SurfElements);

% Step 3: Solve PDE
D = 1; alpha = 1; bcond = 'neu';
% Neumann
f = @(P) 8*(1-2*(P(:,1).^2 + P(:,2).^2)) + (1- (P(:,1).^2 + P(:,2).^2)).^2;
% Dirichlet 
% f = @(P) 5 - P(:,1).^2 - P(:,2).^2;
u = solver_elliptic_bulk(D, alpha, f, P, M, K, R, bcond);

% STEP 4: Post-processing
figure, set(gcf,'color','white')
plot_bulk2d(BulkElements, u, '$u$');
% Neumann
u_exact = (1- (P(:,1).^2 + P(:,2).^2)).^2;
% Dirichlet
% u_exact = 1- (P(:,1).^2 + P(:,2).^2);
L2_relative_err = compute_error(C,[],u,u_exact,[],[]);