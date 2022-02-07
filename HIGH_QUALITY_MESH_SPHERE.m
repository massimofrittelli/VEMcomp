close all
clearvars

% GENERATE 3D MESH WITH DISTMESH
% fd = @(p) sqrt(sum(p.^2,2))-1;
% [p,t] = distmeshnd(fd, @huniform, 0.1, [-1,-1,-1;1,1,1],[]);

load('sphere_tri_4385.mat');

[K,M,KS,MS,ES,R] = matricesTet(PB,EB);

P = PB;
EGamma = ES;

figure
scatter3(PB(:,1), PB(:,2), PB(:,3), 5, 'filled')

figure
trisurf(ES, PB(:,1), PB(:,2), PB(:,3), 'FaceColor', 'none')