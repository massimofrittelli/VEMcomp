%% VEMCOMP:  Example 6.2.1 surface-only linear PARABOLIC PDE eq. (32)
% Linear parabolic problem on a spherical surface with exact solution 
clear all
disp('Linear parabolic problem on a spherical surface')

% STEP 1: Generate mesh
level_fun = @(P) P(:,1).^2 + P(:,2).^2 + P(:,3).^2 -1;
range = [-1,1; -1,1; -1,1];
Nx = 30; tol = 1e-6; xcut = -0.3;
[P,h,BulkElements,SurfElements,ElementsPlot] = ...
    generate_mesh3d(level_fun,range,Nx,tol,xcut);

% STEP 2: Matrix assembly
[K,M,C,KS,MS,CS,R] = assembly3d(P,BulkElements,SurfElements);

% STEP 3: Solve PDE
g = {@(u,P,t) 13*P(:,1).*P(:,2).*P(:,3)*exp(t)};
v0 = R'*(P(:,1).*P(:,2).*P(:,3));
D = 1; T = 1; tau = 1e-4;
v = solver_parabolic_surf(D, g, P, MS, KS, R, T, tau, v0);

% STEP 4: Post-processing
figure, set(gcf,'Color','white')
plot_surf3d(P,R,SurfElements,v,'Linear parabolic problem $v$')
xlim([-0.5,1]) % to cut the surface
v_exact = R'*(P(:,1).*P(:,2).*P(:,3))*exp(1);
L2_relative_err = compute_error([],MS,[],[],v,v_exact);

disp('L2 relative error'),
disp(L2_relative_err)