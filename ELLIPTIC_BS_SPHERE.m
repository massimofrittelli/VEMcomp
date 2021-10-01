% Script to solve an elliptic B-S problem on the sphere
alpha = 1;
beta = 2;

load('sphere41.mat')
N = length(P); % Overall amount of nodes
NGamma = length(MS); % Amount of boundary nodes

esol_u = @(P) P(:,1).*P(:,2).*P(:,3);
esol_v = @(P) 2*P(:,1).*P(:,2).*P(:,3);
f_u = @(P) P(:,1).*P(:,2).*P(:,3);
f_v = @(P) 29*P(:,1).*P(:,2).*P(:,3);

RM = spalloc(N, NGamma, NGamma);
RM(boundarynode, 1:NGamma) = speye(NGamma); % reduction matrix

tic
numsol = [K+M+alpha*RM*MS*RM', -beta*RM*MS; -alpha*MS*RM', KS+(beta+1)*MS]\[M*f_u(P); MS*f_v(P(boundarynode,:))];
u = numsol(1:N,1);
v = numsol(N+1:end,1);
toc
es_u = esol_u(P);
es_v = esol_v(P(boundarynode,:));
err_u = u - es_u;
err_v = v - es_v;
L2err_u = sqrt(err_u'*M*err_u);
L2err_v = sqrt(err_v'*MS*err_v);
L2err_product = sqrt(err_u'*M*err_u + err_v'*MS*err_v);

% normsol = sqrt(es_u'*M*es_u);
% L2errel = L2err/normsol;

% Plotting Numerical Solution - Bulk Component u
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

% Plotting Numerical Solution - Surface Component v
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
