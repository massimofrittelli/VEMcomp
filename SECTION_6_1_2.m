%% VEMCOMP:  Example 6.1.2 bulk-only ELLIPTIC PDE eq. (31)
% Poisson problem in 3D on a spherical domain with zero
% Dirichlet boundary conditions and exact solution 
clear all
disp('3D Poisson problem on a spherical domain, zero Dirichlet bc')

% STEP 1: Generate mesh
level_fun = @(P) P(:,1).^2 + P(:,2).^2 + P(:,3).^2 -1;
range = [-1,1; -1,1; -1,1];
Nx = 30; tol = 1e-6; xcut = -0.3;
[P,h,BulkElements,SurfElements,ElementsPlot] = ...
    generate_mesh3d(level_fun,range,Nx,tol,xcut);

% STEP 2: Matrix assembly
[K,M,C,KS,MS,CS,R] = assembly3d(P,BulkElements,SurfElements);

% STEP 3: Solve PDE
D = 1;   alpha = 1; bcond = 'dir';
f = @(P) 7 - (P(:,1).^2 + P(:,2).^2 + P(:,3).^2);
u = solver_elliptic_bulk(D, alpha, f, P, M, K, R, bcond);

% STEP 4: Post-processing
figure, set(gcf,'color','white')
plot_bulk3d(ElementsPlot, u, 'Poisson problem $u$')
colormap cool % Default colormap is parula
u_exact = 1- (P(:,1).^2 + P(:,2).^2 + P(:,3).^2);
L2_relative_err = compute_error(C,[],u,u_exact,[],[]);

disp('L2 relative error'),
disp(L2_relative_err)