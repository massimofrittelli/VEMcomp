% DESCRIPTION - Solves the Laplace equation on the 3D sphere with NON
% homogeneous Neumann boundary conditions with VEM on a polyhefral mesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clearvars

load('sphere41.mat')
N = length(PB);

NGamma = length(MS); % Amount of boundary nodes
R = spalloc(N, NGamma, NGamma);
R(boundarynode, 1:NGamma) = speye(NGamma); % reduction matrix

esol = @(P) P(:,1).*P(:,2).*P(:,3);
f = @(P)    P(:,1).*P(:,2).*P(:,3);
r = @(P)  3*P(:,1).*P(:,2).*P(:,3);
tic
u = (K+M)\(M*f(PB) + R*MS*r(R'*PB));
toc
es = esol(PB);
err = u - es;
L2err = sqrt(err'*M*err);
normsol = sqrt(es'*M*es);
L2errel = L2err/normsol;

% Plotting Exact solution
figure
indsol = PB(:,1) >= -0.5;
chull = convhull(PB(indsol,1), PB(indsol,2), PB(indsol,3));
set(gcf,'Color','white')
trisurf(chull, PB(indsol,1), PB(indsol,2), PB(indsol,3), es(indsol,1),'EdgeColor', 'none', 'FaceColor', 'interp')
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
indsol = PB(:,1) >= -0.5;
chull = convhull(PB(indsol,1), PB(indsol,2), PB(indsol,3));
set(gcf,'Color','white')
trisurf(chull, PB(indsol,1), PB(indsol,2), PB(indsol,3), u(indsol,1),'EdgeColor', 'none', 'FaceColor', 'interp')
view(3)
set(gca,'FontSize',18)
xlabel('x')
ylabel('y')
zlabel('z','rot',0)
axis equal
title('Numerical solution')
colorbar
