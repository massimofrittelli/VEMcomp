% DESCRIPTION: Solves BULK_ONLY linear heat equation on the 3D sphere, with
% non-homogeneous Neumann boundary conditions, plots solutions and computes
% errors.

close all
clearvars

T = 1;
tau = .01;

% Generating mesh
fun = @(P) P(:,1).^2 + P(:,2).^2 + P(:,3).^2 - 1;
range = [-1,1; -1,1; -1,1];
tol = 1e-3;
Nx = 10;
xcut = -0.5;
[P, h, BulkElements, SurfElements, ElementsPlot] = generate_mesh3d(fun, range, Nx, tol, xcut);

% Assembling matrices
[K, M, C, KS, MS, CS, R] = assembly3d(P, BulkElements, SurfElements);

N = length(P);

NGamma = length(MS); % Amount of boundary nodes

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

% Plotting numerical solution
figure
set(gcf,'color','white')
for i=1:length(ElementsPlot)
   plot(ElementsPlot(i), u(ElementsPlot(i).Pind)); 
   hold on
end
colormap jet
axis equal tight
xlabel('x')
ylabel('y','rot',0)
zlabel('z','rot',0)
title('u')
set(gca,'FontSize',18)
colorbar

% Plotting exact solution
figure
set(gcf,'color','white')
for i=1:length(ElementsPlot)
   plot(ElementsPlot(i), esol(P(ElementsPlot(i).Pind, :), T)); 
   hold on
end
colormap jet
axis equal tight
xlabel('x')
ylabel('y','rot',0)
zlabel('z','rot',0)
set(gca,'FontSize',18)
title('u_{exact}')
colorbar