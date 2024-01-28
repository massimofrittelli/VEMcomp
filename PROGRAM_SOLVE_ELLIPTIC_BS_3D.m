% DESCRIPTION - Solves an elliptic B-S toy model on the sphere, plots solution and
% computes errors.

close all
clearvars

% STEP 1: Generate mesh
level_fun = @(P) P(:,1).^2 + P(:,2).^2 + P(:,3).^2 - 1;
range = [-1,1; -1,1; -1,1];
Nx = 10; tol = 1e-6; xcut = -0.3;
[P,h,BulkElements,SurfElements,ElementsPlot] = ...
    generate_mesh3d(level_fun, range, Nx, tol, xcut);

% STEP 2: Matrix assembly
[K,M,C,KS,MS,CS,R] = assembly3d(P,BulkElements,SurfElements);

% STEP 3: Solve PDE
dOmega = 1; dGamma = 1;
alpha = 1; beta = 2; gamma = -1; delta = 2;
% right-hand-side of 1st eq. in (35)
f = @(P) P(:,1).*P(:,2).*P(:,3);
% right-hand-side of 2nd eq. in (35)
g = @(P) 29*P(:,1).*P(:,2).*P(:,3);
[u,v] = solver_elliptic_bulk_surf(dOmega, dGamma, alpha, beta, f, g, P, M, MS, K, KS, R, gamma, delta);

% STEP 4: Post-processing
figure, set(gcf, 'Color','white')
subplot(121) % plot bulk component u
plot_bulk3d(ElementsPlot, u, '$u$')
subplot(122) % plot surf component v
plot_surf3d(P,R,SurfElements,v,'$v$')
xlim([-0.5,1])
% hold on
% Ccirc = [-0.5, 0, 0];   % Center of circle 
% Rcirc = sqrt(3)/2;      % Radius of circle 
% teta = 0:0.01:2*pi;
% xcirc = Ccirc(1) + zeros(size(teta)) ;
% ycirc = Ccirc(2) + Rcirc*cos(teta);
% zcirc = Ccirc(3) + Rcirc*sin(teta) ;
% plot3(xcirc,ycirc,zcirc,'k','LineWidth',1.5)
u_exact = P(:,1).*P(:,2).*P(:,3);
v_exact = 2*R'*(P(:,1).*P(:,2).*P(:,3));
L2_err_rel = compute_error(C,MS,u,u_exact,v,v_exact);


% TOY PROBLEM 2 - 2 EIGENMODES
% esol_u = @(P) P(:,1).*P(:,2).*P(:,3) - P(:,1).*P(:,2);
% esol_v = @(P) 2*P(:,1).*P(:,2).*P(:,3) - 3/2*P(:,1).*P(:,2);
% f_u = @(P) P(:,1).*P(:,2).*P(:,3) - P(:,1).*P(:,2);
% f_v = @(P) 29*P(:,1).*P(:,2).*P(:,3) - 25/2*P(:,1).*P(:,2);