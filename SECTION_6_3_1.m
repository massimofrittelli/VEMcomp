%% VEMCOMP: Section 6.3 bulk-surface problems in 3D
%     Example 6.3.1 Elliptic bulk-surface problem on the sphere with exact solution
clear all
disp('Elliptic bulk-surface problem on the sphere')

%5 STEP 1: Generate mesh
level_fun = @(P) P(:,1).^2 + P(:,2).^2 + P(:,3).^2 - 1;
range = [-1,1; -1,1; -1,1];
Nx = 20; 
tol = 1e-6; xcut = -0.3;
[P,h,BulkElements,SurfElements,ElementsPlot] = ...
    generate_mesh3d(level_fun, range, Nx, tol, xcut);

%% STEP 2: Matrix assembly
[K,M,C,KS,MS,CS,R] = assembly3d(P,BulkElements,SurfElements);

%% STEP 3: Solve PDE
dOmega = 1; dGamma = 1;
alpha = 1; beta = 2; gamma = -1; delta = 2;
% right-hand-side of 1st eq. in (35)
f = @(P) P(:,1).*P(:,2).*P(:,3);
% right-hand-side of 2nd eq. in (35)
g = @(P) 29*P(:,1).*P(:,2).*P(:,3);
[u,v] = solver_elliptic_bulk_surf(dOmega, dGamma, alpha, beta, f, g, P, M, MS, K, KS, R, gamma, delta);

%% STEP 4: Post-processing
figure, set(gcf, 'Color','white')
subplot(121) % plot bulk component u
plot_bulk3d(ElementsPlot, u, 'Bulk solution $u$')
subplot(122) % plot surf component v
plot_surf3d(P,R,SurfElements,v,'Surface solution $v$')
xlim([-0.5,1])
u_exact = P(:,1).*P(:,2).*P(:,3);
v_exact = 2*R'*(P(:,1).*P(:,2).*P(:,3));
L2_relative_err = compute_error(C,MS,u,u_exact,v,v_exact);

disp('L2 relative error'),
disp(L2_relative_err)