% Script to test a parabolic coupled  bulk-surface problem on the 3D CUBE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% u_t - d_u*Laplace(u) = au,                                   x in \Omega
% v_t - d_v*LaplaceBeltrami(v) + \nabla u \cdot \nu = bu+cv,   x in \Gamma
% \nabla u \cdot \nu = du+ev,                                  x in \Gamma
%
% \Omega = [-1/2, 1/2]^3
% \Gamma = \partial \Omega
%
% Exact solution:
%
% u(x,y,z,t) = cos(pi x)cos(pi y)cos(pi z) exp(-alpha t),  (x,y,z) \in \Omega
% v(x,y,z,t) = cos(pi x)cos(pi y)cos(pi z) exp(-alpha t),  (x,y,z) \in \Gamma
%
% if the model parameters fulfil
%
% a = - alpha + 3 d_u pi^2
% b     arbitrary
% c = - alpha + 2 d_v pi^2 - pi
% d     arbitrary
% e = - pi
%
% Note: to test
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


close all
clearvars

% Set solution parameters
alpha = 1;

% Set free model parameters
du = 1;
dv = 1;
b = 2;
d = 1;

% Dependent model parameters (DO NOT TOUCH)
a = - alpha + 3*du*pi^2;
c = - alpha + 2*dv*pi^2 - pi;
e = - pi;

% Set time discretisation
T = 1;
%tau = 1.5625e-02;
tau = 1e-3;

% Set space discretisation
load('cube17.mat')
P = P - 0.5;

% END OF INPUT

N = length(P);
NGamma = length(MGamma); % Amount of boundary nodes

esol_u = @(P,t) cos(pi*P(:,1)).*cos(pi*P(:,2)).*cos(pi*P(:,3))*exp(-alpha*t);
esol_v = @(P,t) cos(pi*P(:,1)).*cos(pi*P(:,2)).*cos(pi*P(:,3))*exp(-alpha*t);

tic
NT = ceil(T/tau);
u =  esol_u(P,0);
v =  esol_v(A'*P,0);

%for i=0:NT-1
for i=0
   
   % IMEX FULLY AUTONOMOUS - shows restrictive conditions on timestep
   unext = (M+tau*du*K)\(M*(1+tau*a)*u + tau*du*A*MGamma*(d*A'*u+e*v)); % - not working
   vnext = (MGamma+tau*dv*KGamma)\(MGamma*(v+tau*((b-d)*A'*u + (c-e)*v)));
   
   u = unext;
   v = vnext;

end
toc

es_u = esol_u(P,T);
es_v = esol_v(A'*P,T);
err_u = u - es_u;
err_v = v - es_v;
L2err_u = sqrt(err_u'*M*err_u);
L2err_v = sqrt(err_v'*MGamma*err_v);
L2err_prod = sqrt(err_u'*M*err_u + err_v'*MGamma*err_v);

% Plotting Exact solution in the bulk
figure
indsol = P(:,1) >= 0;
chull = convhull(P(indsol,1), P(indsol,2), P(indsol,3));
set(gcf,'Color','white')
trisurf(chull, P(indsol,1), P(indsol,2), P(indsol,3), es_u(indsol,1),'EdgeColor', 'none', 'FaceColor', 'interp')
view(3)
set(gca,'FontSize',18)
xlabel('x')
ylabel('y')
zlabel('z','rot',0)
axis equal
title('Exact solution')
colorbar

% Plotting Numerical solution in the bulk
figure
indsol = P(:,1) >= 0;
chull = convhull(P(indsol,1), P(indsol,2), P(indsol,3));
set(gcf,'Color','white')
trisurf(chull, P(indsol,1), P(indsol,2), P(indsol,3), u(indsol,1),'EdgeColor', 'none', 'FaceColor', 'interp')
view(3)
set(gca,'FontSize',18)
xlabel('x')
ylabel('y')
zlabel('z','rot',0)
axis equal
title('Bulk solution u')
colorbar

% Plotting Numerical Solution on the surface
figure
set(gcf,'Color','white')
trisurf(Egamma, P(:,1), P(:,2), P(:,3), [zeros(N-NGamma,1); v], 'EdgeColor', 'none', 'FaceColor', 'interp')
view(3)
set(gca,'FontSize',18)
xlabel('x')
ylabel('y')
zlabel('z','rot',0)
axis equal
xlim([-0.5,1])
title('Surface solution v')
colorbar