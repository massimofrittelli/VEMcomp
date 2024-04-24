% DESCRIPTION - Solves an elliptic B-S toy model on the sphere, plots solution and
% computes errors.

close all
clearvars

% Generating mesh
fun = @(P) P(:,1).^2 + P(:,2).^2 + P(:,3).^2 - 1;
range = [-1,1; -1,1; -1,1];
tol = 1e-3;
Nx = 10;
xcut = -0.5;
[P, h, BulkElements, SurfElements, ElementsPlot] = generate_mesh3d(fun, range, Nx, tol, xcut);

% Assembling matrices
[K, M, C, KS, MS, CS, R] = assembly3d(P, BulkElements, SurfElements);

T = 1;
tau = 1e-2;

[u, v, L2err_prod] = solver_parabolic3d_bs_sphere(SurfElements, ElementsPlot, P, M, K, MS, KS, R, T, tau);