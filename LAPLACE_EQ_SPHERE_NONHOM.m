% Script to test the Laplace equation on the 3D sphere, with
% NON-HOMOGENEOUS Neumann boundary conditions

close all
clearvars

load('sphere41.mat')
N = length(P);

NGamma = length(MS); % Amount of boundary nodes
R = spalloc(N, NGamma, NGamma);
R(boundarynode, 1:NGamma) = speye(NGamma); % reduction matrix

esol = @(P) P(:,1).*P(:,2).*P(:,3);
f = @(P)    P(:,1).*P(:,2).*P(:,3);
r = @(P)  3*P(:,1).*P(:,2).*P(:,3);
tic
u = (K+M)\(M*f(P) + R*MS*r(R'*P));
toc
es = esol(P);
err = u - es;
L2err = sqrt(err'*M*err);
normsol = sqrt(es'*M*es);
L2errel = L2err/normsol;

% Plotting Exact solution
figure
indsol = P(:,1) >= -0.5;
chull = convhull(P(indsol,1), P(indsol,2), P(indsol,3));
set(gcf,'Color','white')
trisurf(chull, P(indsol,1), P(indsol,2), P(indsol,3), es(indsol,1),'EdgeColor', 'none', 'FaceColor', 'interp')
view(3)
set(gca,'FontSize',18)
xlabel('x')
ylabel('y')
zlabel('z','rot',0)
axis equal
title('Exact solution')
colorbar

% Plotting Numerical solution
figure
indsol = P(:,1) >= -0.5;
chull = convhull(P(indsol,1), P(indsol,2), P(indsol,3));
set(gcf,'Color','white')
trisurf(chull, P(indsol,1), P(indsol,2), P(indsol,3), u(indsol,1),'EdgeColor', 'none', 'FaceColor', 'interp')
view(3)
set(gca,'FontSize',18)
xlabel('x')
ylabel('y')
zlabel('z','rot',0)
axis equal
title('Numerical solution')
colorbar
