% Generating mesh
fun = @(P) P(:,1).^2 + P(:,2).^2 + P(:,3).^2 - 1;
range = [-1,1; -1,1; -1,1];
tol = 1e-3;
Nx = 10;
xcut = -0.5;
[P, h, BulkElements, SurfaceElements, ElementsPlot] = generate_mesh_3d(fun, range, Nx, tol, xcut);

% Assembling matrices
[K, M, C, KS, MS, CS, R] = assembly_3d(P, BulkElements, SurfaceElements);

% Solving PDE
f = {@(P) 13*P(:,1).*P(:,2).*P(:,3)};
D = 1;
alpha = 1;
v = solver_elliptic_surf(1, D, alpha, f, P, MS, KS, R);

% Compute error
v_exact = R'*(P(:,1).*P(:,2).*P(:,3));
L2_relative_err = compute_error([],MS,[],[],v,v_exact);

% Plotting numerical solution
figure
set(gcf,'color','white')
plot_surf_3d(P,R,SurfaceElements,v,'$v$')
xlim([xcut, 1])

% Plotting exact solution
esol = @(P) P(:,1).*P(:,2).*P(:,3);
figure
set(gcf,'color','white')
plot_surf_3d(P,R,SurfaceElements,esol(R'*P),'$v_{exact}$')
xlim([xcut, 1])