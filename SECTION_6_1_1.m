%% VEMCOMP:  Example 6.1.1 bulk-only ELLIPTIC PDE eq. (30)
% Poisson problem in 2D on a circular domain with zero
% Neumann boundary conditions and exact solution 
clear all
disp('2D Poisson problem on a circular domain, zero Neumann bc')

%% STEP 1: Generate mesh
level_fun = @(P) P(:,1).^2 + P(:,2).^2 -1;
range = [-1,1; -1,1]; Nx = 40; tol = 1e-6;
[P, h, BulkElements, SurfElements] = generate_mesh2d(level_fun, range, Nx, tol);

%% STEP 2: Matrix assembly
[K,C,M,KS,MS,R] = assembly2d(P, BulkElements, SurfElements);

%% STEP3: Solve PDE
D = 1; alpha = 1; bcond = 'neu';
f = @(P) 8*(1-2*(P(:,1).^2 + P(:,2).^2)) + (1- (P(:,1).^2 + P(:,2).^2)).^2;
u = solver_elliptic_bulk(D, alpha, f, P, M, K, R, bcond);

%% STEP 4: Post-processing
figure, set(gcf,'color','white')
plot_bulk2d(BulkElements, u, 'Poisson Problem $u$');
u_exact = (1- (P(:,1).^2 + P(:,2).^2)).^2;
L2_relative_err = compute_error(C,[],u,u_exact,[],[]);

disp('L2 relative error'),
disp(L2_relative_err)