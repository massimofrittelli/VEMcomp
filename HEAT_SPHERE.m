% DESCRIPTION: Solves BULK_ONLY linear heat equation on the 3D sphere, with
% non-homogeneous Neumann boundary conditions, plots solutions and computes
% errors.

close all
clearvars

T = 1;
tau = .01;

load('sphere21.mat')
N = length(PB);

NGamma = length(MS); % Amount of boundary nodes
R = spalloc(N, NGamma, NGamma);
R(boundarynode, 1:NGamma) = speye(NGamma); % reduction matrix

esol = @(P,t) P(:,1).*P(:,2).*P(:,3)*exp(t);
f = @(u)      u;
r = @(P,t)    3*P(:,1).*P(:,2).*P(:,3)*exp(t);

tic
NT = ceil(T/tau);
u = esol(PB,0);
for i=0:NT-1
   u = (M+tau*K)\(M*(u+tau*f(u)) + tau*R*MS*r(R'*PB,i*tau)); 
end
toc

es = esol(PB,T);
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