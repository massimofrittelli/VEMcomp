% Script to test a parabolic coupled  bulk.surface problem on the 3D sphere

close all
clearvars

a = 1;
b = 2;
du = 1;
dv = 1;
T = 1;
tau = 1.5625e-02;

load('sphere41.mat')
N = length(P);

NGamma = length(MS); % Amount of boundary nodes
R = spalloc(N, NGamma, NGamma);
R(boundarynode, 1:NGamma) = speye(NGamma); % reduction matrix

esol_u = @(P,t) a*P(:,1).*P(:,2).*P(:,3)*exp(t);
esol_v = @(P,t) b*P(:,1).*P(:,2).*P(:,3)*exp(t);
q = @(P,t) a*P(:,1).*P(:,2).*P(:,3)*exp(t);
s = @(P,t) 3*a*P(:,1).*P(:,2).*P(:,3)*exp(t);
r_minus_s = @(P,t) (12*dv+1)*b*P(:,1).*P(:,2).*P(:,3)*exp(t);

tic
NT = ceil(T/tau);
u = esol_u(P,0);
v = esol_v(R'*P,0);
for i=0:NT-1
   unext = (M+tau*du*K)\(M*(u+tau*q(P,i*tau)) + tau*du*R*MS*s(R'*P,i*tau));
   vnext = (MS+tau*dv*KS)\(MS*(v+tau*r_minus_s(R'*P,i*tau)));
   u = unext;
   v = vnext;
end
toc

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