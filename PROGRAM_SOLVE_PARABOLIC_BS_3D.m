% DESCRIPTION - Solves an elliptic B-S toy model on the sphere, plots solution and
% computes errors.

close all

%load('mesh_sphere_marchcub_Nx31.mat')

T = 1;
tau = 1e-2;

[u, v, L2err_prod] = solver_parabolic3d_bs_sphere(SurfaceElements, ElementsPlot, P, M, K, MS, KS, R, T, tau);