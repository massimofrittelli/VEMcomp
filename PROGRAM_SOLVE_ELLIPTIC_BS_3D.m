% DESCRIPTION - Solves an elliptic B-S toy model on the sphere, plots solution and
% computes errors.

close all
%clearvars

% Script to solve an elliptic B-S problem on the sphere
alpha = 1;
beta = 2;

%load('mesh_sphere_marchcub_Nx30.mat')
%load('mesh_sphere31.mat')
%load('mesh_sphere30_tol1e-6.mat')
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
u_exact = esol_u(P);
v_exact = esol_v(R'*P);
err_u = u - u_exact;
err_v = v - v_exact;

% L2err_product_abs = sqrt(err_u'*M*err_u + err_v'*MS*err_v);
% H1err_product_abs = sqrt(err_u'*(K+M)*err_u + err_v'*(KS+MS)*err_v);
% L2norm_product = sqrt(u_exact'*M*u_exact + v_exact'*MS*v_exact);
% H1norm_product = sqrt(u_exact'*(K+M)*u_exact + v_exact'*(KS+MS)*v_exact);
% L2err_product_rel = L2err_product_abs/L2norm_product;
% H1err_product_rel = H1err_product_abs/H1norm_product;

L2_err_rel = compute_error(C,MS,u,u_exact,v,v_exact);


%xcut = range(1,1) + h/sqrt(3)*2;

% Plotting Numerical Solution
figure
set(gcf, 'Color','white')
% Bulk Component u
subplot(121)
plot_bulk_3d(ElementsPlot, u, '$u$')
% Surface Component v
subplot(122)
plot_surf_3d(P,R,SurfElements,v,'$v$')
xlim([xcut,1])
hold on
Ccirc = [-0.5, 0, 0];   % Center of circle 
Rcirc = sqrt(3)/2;      % Radius of circle 
teta = 0:0.01:2*pi;
xcirc = Ccirc(1) + zeros(size(teta)) ;
ycirc = Ccirc(2) + Rcirc*cos(teta);
zcirc = Ccirc(3) + Rcirc*sin(teta) ;
plot3(xcirc,ycirc,zcirc,'k','LineWidth',1.5)