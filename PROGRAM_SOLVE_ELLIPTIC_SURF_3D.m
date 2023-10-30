% Generating mesh
fun = @(P) P(:,1).^2 + P(:,2).^2 + P(:,3).^2 - 1;
range = [-1,1; -1,1; -1,1];
tol = 1e-3;
Nx = 20;
[P, h, BulkElements, SurfaceElements, ElementsPlot] = generate_mesh_3d(fun, range, Nx, tol, -0.5);

% Assembling matrices
[K, M, C, KS, MS, CS, R] = assembly_3d(P, BulkElements, SurfaceElements);

% Solving PDE
f = {@(P) 13*P(:,1).*P(:,2).*P(:,3)};
D = 1;
alpha = 1;
u = solver_elliptic_surf(1, D, alpha, f, P, MS, KS, R);

% Plotting numerical solution
figure
set(gcf,'color','white')
trisurf(SurfaceElements, P(:,1), P(:,2), P(:,3), R*u, 'FaceColor', 'interp');
colormap jet
axis equal tight
xlabel('x')
ylabel('y','rot',0)
zlabel('z','rot',0)
title('u_h')
set(gca,'FontSize',18)
colorbar

esol = @(P) P(:,1).*P(:,2).*P(:,3);
figure
set(gcf,'color','white')
trisurf(SurfaceElements, P(:,1), P(:,2), P(:,3), esol(P), 'FaceColor', 'interp');
colormap jet
axis equal tight
xlabel('x')
ylabel('y','rot',0)
zlabel('z','rot',0)
set(gca,'FontSize',18)
title('u_{exact}')
colorbar