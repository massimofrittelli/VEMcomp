% DESCRIPTION - Solves an elliptic B-S toy model on the sphere, plots solution and
% computes errors.

close all

% Script to solve an elliptic B-S problem on the sphere
alpha = 1;
beta = 2;

load('mesh_sphere_marchcub_Nx30.mat')
N = length(P); % Overall amount of nodes
NGamma = length(MS); % Amount of boundary nodes

% TOY PROBLEM 1 - 1 EIGENMODE
esol_u = @(P) P(:,1).*P(:,2).*P(:,3);
esol_v = @(P) 2*P(:,1).*P(:,2).*P(:,3);
f_u = @(P) P(:,1).*P(:,2).*P(:,3);
f_v = @(P) 29*P(:,1).*P(:,2).*P(:,3);

% TOY PROBLEM 2 - 2 EIGENMODES
% esol_u = @(P) P(:,1).*P(:,2).*P(:,3) - P(:,1).*P(:,2);
% esol_v = @(P) 2*P(:,1).*P(:,2).*P(:,3) - 3/2*P(:,1).*P(:,2);
% f_u = @(P) P(:,1).*P(:,2).*P(:,3) - P(:,1).*P(:,2);
% f_v = @(P) 29*P(:,1).*P(:,2).*P(:,3) - 25/2*P(:,1).*P(:,2);

tic
numsol = [K+M+alpha*R*MS*R', -beta*R*MS; -alpha*MS*R', KS+(beta+1)*MS]\[M*f_u(P); MS*f_v(R'*P)];
u = numsol(1:N,1);
v = numsol(N+1:end,1);
toc
es_u = esol_u(P);
es_v = esol_v(R'*P);
err_u = u - es_u;
err_v = v - es_v;
% L2err_u = sqrt(err_u'*M*err_u);
% L2err_v = sqrt(err_v'*MS*err_v);
% H1err_u = sqrt(err_u'*K*err_u);
% H1err_v = sqrt(err_v'*KS*err_v);

L2err_product_abs = sqrt(err_u'*M*err_u + err_v'*MS*err_v);
H1err_product_abs = sqrt(err_u'*(K+M)*err_u + err_v'*(KS+MS)*err_v);
L2norm_product = sqrt(es_u'*M*es_u + es_v'*MS*es_v);
H1norm_product = sqrt(es_u'*(K+M)*es_u + es_v'*(KS+MS)*es_v);
L2err_product_rel = L2err_product_abs/L2norm_product;
H1err_product_rel = H1err_product_abs/H1norm_product;


xcut = range(1,1) + h/sqrt(3)*2;


% Plotting Numerical Solution - Bulk Component u
figure
set(gcf, 'Color','white')

subplot(121)
hold on
for i=1:length(ElementsPlot)
       plot(ElementsPlot(i), u(ElementsPlot(i).Pind)); 
end
view(3)
set(gca,'FontSize',18)
xlabel('x')
ylabel('y')
zlabel('z','rot',0)
axis equal
title('u')
colorbar
colormap parula

% Plotting Numerical Solution - Surface Component v
subplot(122)
trisurf(SurfaceElements, P(:,1), P(:,2), P(:,3), R*v, 'FaceColor', 'interp')
view(3)
set(gca,'FontSize',18)
xlabel('x')
ylabel('y')
zlabel('z','rot',0)
axis equal
xlim([-0.5,1])
title('v')
colorbar
colormap parula
hold on
Ccirc = [-0.5, 0, 0];   % Center of circle 
Rcirc = sqrt(3)/2;      % Radius of circle 
teta = 0:0.01:2*pi;
xcirc = Ccirc(1) + zeros(size(teta)) ;
ycirc = Ccirc(2) + Rcirc*cos(teta);
zcirc = Ccirc(3) + Rcirc*sin(teta) ;
plot3(xcirc,ycirc,zcirc,'k','LineWidth',1.5)