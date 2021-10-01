% Script to test the linear heat equation on the 3D sphere, with
% NON-HOMOGENEOUS Neumann boundary conditions

close all
clearvars

T = 1;
tau = .01;

load('sphere21.mat')
N = length(P);

NGamma = length(MS); % Amount of boundary nodes
R = spalloc(N, NGamma, NGamma);
R(boundarynode, 1:NGamma) = speye(NGamma); % reduction matrix

esol = @(P,t) P(:,1).*P(:,2).*P(:,3)*exp(t);
f = @(u)      u;
r = @(P,t)    3*P(:,1).*P(:,2).*P(:,3)*exp(t);

tic
NT = ceil(T/tau);
u = esol(P,0);
for i=0:NT-1
   u = (M+tau*K)\(M*(u+tau*f(u)) + tau*R*MS*r(R'*P,i*tau)); 
end
toc

es = esol(P,T);
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