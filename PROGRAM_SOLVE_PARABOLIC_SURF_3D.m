% DESCRIPTION - Solves a parabolic surface toy model on the sphere, plots
% solution and computes errors.

close all

load('mesh_sphere_marchcub_Nx31.mat')

T = 1;
tau = 1e-4;

[v, L2err] = solver_parabolic3d_surf_sphere(SurfaceElements, R, P, MS, KS, T, tau);