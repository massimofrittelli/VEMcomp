% Script to test a parabolic coupled  bulk-surface problem on the 3D sphere

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% u_t - d_u*Laplace(u) = au,                                   x in \Omega
% v_t - d_v*LaplaceBeltrami(v) + \nabla u \cdot \nu = bu+cv,   x in \Gamma
% \nabla u \cdot \nu = du+ev,                                  x in \Gamma
%
% \Omega = unit sphere
%
% Exact solution:
%
% u(x,t) =      xyz exp(-alpha t)
% v(x,t) = beta xyz exp(-alpha t)
%
% if the model parameters fulfil
%
% a = - alpha
% b = 3 + beta(12d_v -alpha-c)
% d = 3 - beta e
%
% Note: if the solution is exponentially-in-time increasing instead of
% decreasing, it causes stability issues with the method. 
% The solution blows up. The stability condition probably becomes too
% restrictive.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


close all
clearvars

% Set solution parameters
alpha = 1;
beta = 2;

% Set free model parameters
du = 1;
dv = .5;
c = 2;
e = 1;

% Dependent model parameters (DO NOT TOUCH)
a = - alpha;
b = 3 + beta*(12*dv -alpha -c);
d = 3 - beta*e;

% Set time discretisation
T = 1;
tau = 1.5625e-02;
%tau = 1e-3;

% Set space discretisation
load('sphere21.mat')

% END OF INPUT

N = length(P);
NGamma = length(MS); % Amount of boundary nodes
R = spalloc(N, NGamma, NGamma);
R(boundarynode, 1:NGamma) = speye(NGamma); % reduction matrix

esol_u = @(P,t) P(:,1).*P(:,2).*P(:,3)*exp(-alpha*t);
esol_v = @(P,t) beta*P(:,1).*P(:,2).*P(:,3)*exp(-alpha*t);

tic
NT = ceil(T/tau);
u = esol_u(P,0);
v = esol_v(R'*P,0);
for i=0:NT-1
   %unext = (M+tau*du*(K-(1+beta)*R*MS*R'))\(M*(1+tau*a)*u); % V3 - not working
   unext = (M+tau*du*K)\(M*(1+tau*a)*u + tau*du*R*MS*3*esol_u(R'*P,i*tau)); % V2 - working
   %unext = (M+tau*du*K)\(M*(1+tau*a)*u + tau*du*R*MS*(d*R'*u+e*v)); % V1 - not working
   vnext = (MS+tau*dv*KS)\(MS*(v+tau*((b-d)*R'*u + (c-e)*v)));
   u = unext;
   v = vnext;
end
toc

% TODO: try and eliminate the stability part from the mass matrix on the
% RHSs

es_u = esol_u(P,T);
es_v = esol_v(R'*P,T);
err_u = u - es_u;
err_v = v - es_v;
L2err_u = sqrt(err_u'*M*err_u);
L2err_v = sqrt(err_v'*MS*err_v);
L2err_prod = sqrt(err_u'*M*err_u + err_v'*MS*err_v);

% Plotting Exact solution in the bulk
figure
indsol = P(:,1) >= -0.5;
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