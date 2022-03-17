% DESCRIPTION - Solves an elliptic B-S toy model on the sphere, plots solution and
% computes errors.

close all

% Script to solve an elliptic B-S problem on the sphere
alpha = 1;
beta = 2;

load('sphere6.mat')
N = length(PB); % Overall amount of nodes
NGamma = length(MS); % Amount of boundary nodes

% TOY PROBLEM 1 - 1 EIGENMODE
% esol_u = @(P) P(:,1).*P(:,2).*P(:,3);
% esol_v = @(P) 2*P(:,1).*P(:,2).*P(:,3);
% f_u = @(P) P(:,1).*P(:,2).*P(:,3);
% f_v = @(P) 29*P(:,1).*P(:,2).*P(:,3);

% TOY PROBLEM 2 - 2 EIGENMODES
esol_u = @(P) P(:,1).*P(:,2).*P(:,3) - P(:,1).*P(:,2);
esol_v = @(P) 2*P(:,1).*P(:,2).*P(:,3) - 3/2*P(:,1).*P(:,2);
f_u = @(P) P(:,1).*P(:,2).*P(:,3) - P(:,1).*P(:,2);
f_v = @(P) 29*P(:,1).*P(:,2).*P(:,3) - 25/2*P(:,1).*P(:,2);

RM = spalloc(N, NGamma, NGamma);
RM(boundarynode, 1:NGamma) = speye(NGamma); % reduction matrix

tic
numsol = [K+M+alpha*RM*MS*RM', -beta*RM*MS; -alpha*MS*RM', KS+(beta+1)*MS]\[M*f_u(PB); MS*f_v(PB(boundarynode,:))];
u = numsol(1:N,1);
v = numsol(N+1:end,1);
toc
es_u = esol_u(PB);
es_v = esol_v(PB(boundarynode,:));
err_u = u - es_u;
err_v = v - es_v;
L2err_u = sqrt(err_u'*M*err_u);
L2err_v = sqrt(err_v'*MS*err_v);
H1err_u = sqrt(err_u'*K*err_u);
H1err_v = sqrt(err_v'*KS*err_v);

L2err_product = sqrt(err_u'*M*err_u + err_v'*MS*err_v);
H1err_product = sqrt(err_u'*K*err_u + err_v'*KS*err_v);

% normsol = sqrt(es_u'*M*es_u);
% L2errel = L2err/normsol;

% Plotting Numerical Solution - Bulk Component u
figure
indsol = PB(:,1) >= -0.5;
chull = convhull(PB(indsol,1), PB(indsol,2), PB(indsol,3));
set(gcf, 'Color','white')
set(gcf, 'Renderer','painters')
%set(gcf, 'Renderer', 'zbuffer')
trisurf(chull, PB(indsol,1), PB(indsol,2), PB(indsol,3), u(indsol,1),'EdgeColor', 'none', 'FaceColor', 'interp')
view(3)
set(gca,'FontSize',18)
xlabel('x')
ylabel('y')
zlabel('z','rot',0)
axis equal
title('$U$', 'interpreter', 'latex')
colorbar
hold on
Ccirc = [-0.5, 0, 0];   % Center of circle 
Rcirc = sqrt(3)/2;      % Radius of circle 
teta = 0:0.01:2*pi;
xcirc = Ccirc(1) + zeros(size(teta)) ;
ycirc = Ccirc(2) + Rcirc*cos(teta);
zcirc = Ccirc(3) + Rcirc*sin(teta) ;
plot3(xcirc,ycirc,zcirc,'k','LineWidth',1.5)
% lightangle(gca,-65,30)
% lighting gouraud

% Plotting Numerical Solution - Surface Component v
figure
set(gcf,'Color','white')
trisurf(EGamma, PB(:,1), PB(:,2), PB(:,3), [zeros(N-NGamma,1); v], 'EdgeColor', 'none', 'FaceColor', 'interp')
view(3)
set(gca,'FontSize',18)
xlabel('x')
ylabel('y')
zlabel('z','rot',0)
axis equal
xlim([-0.5,1])
title('$V$', 'interpreter', 'latex')
colorbar
hold on
Ccirc = [-0.5, 0, 0];   % Center of circle 
Rcirc = sqrt(3)/2;      % Radius of circle 
teta = 0:0.01:2*pi;
xcirc = Ccirc(1) + zeros(size(teta)) ;
ycirc = Ccirc(2) + Rcirc*cos(teta);
zcirc = Ccirc(3) + Rcirc*sin(teta) ;
plot3(xcirc,ycirc,zcirc,'k','LineWidth',1.5)
% lightangle(gca,-65,30)
% lighting gouraud

% Plotting Numerical Error - Bulk Component u
figure
indsol = PB(:,1) >= -0.5;
chull = convhull(PB(indsol,1), PB(indsol,2), PB(indsol,3));
set(gcf, 'Color','white')
set(gcf, 'Renderer','painters')
%set(gcf, 'Renderer', 'zbuffer')
trisurf(chull, PB(indsol,1), PB(indsol,2), PB(indsol,3), es_u(indsol,1)-u(indsol,1),'EdgeColor', 'none', 'FaceColor', 'interp')
view(3)
set(gca,'FontSize',18)
xlabel('x')
ylabel('y')
zlabel('z','rot',0)
axis equal
caxis([min(es_u-u), max(es_u-u)])
title('$u^{-\ell}-U$', 'interpreter', 'latex')
colorbar
hold on
Ccirc = [-0.5, 0, 0];   % Center of circle 
Rcirc = sqrt(3)/2;      % Radius of circle 
teta = 0:0.01:2*pi;
xcirc = Ccirc(1) + zeros(size(teta)) ;
ycirc = Ccirc(2) + Rcirc*cos(teta);
zcirc = Ccirc(3) + Rcirc*sin(teta) ;
plot3(xcirc,ycirc,zcirc,'k','LineWidth',1.5)
% lightangle(gca,-65,30)
% lighting gouraud

% Plotting Numerical Error - Surface Component v
figure
set(gcf,'Color','white')
trisurf(EGamma, PB(:,1), PB(:,2), PB(:,3), [zeros(N-NGamma,1); es_v-v], 'EdgeColor', 'none', 'FaceColor', 'interp')
view(3)
set(gca,'FontSize',18)
xlabel('x')
ylabel('y')
zlabel('z','rot',0)
axis equal
caxis([min(es_v-v), max(es_v-v)])
xlim([-0.5,1])
title('$v^{-\ell}-V$','interpreter', 'latex')
colorbar
hold on
Ccirc = [-0.5, 0, 0];   % Center of circle 
Rcirc = sqrt(3)/2;      % Radius of circle 
teta = 0:0.01:2*pi;
xcirc = Ccirc(1) + zeros(size(teta)) ;
ycirc = Ccirc(2) + Rcirc*cos(teta);
zcirc = Ccirc(3) + Rcirc*sin(teta) ;
plot3(xcirc,ycirc,zcirc,'k','LineWidth',1.5)
% lightangle(gca,-65,30)
% lighting gouraud
