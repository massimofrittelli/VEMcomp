%% VEMCOMP: Section 6.3 bulk-surface problems in 3D
%  Example 6.3.2 Bulk-surface reaction-diffusion system on the torus (eq. 38)
clear all
disp('Bulk-surface reaction-diffusion system on the torus')



% STEP 1: Assign domain and generate mesh
level_fun = @(P) (sqrt(P(:,1).^2 + P(:,2).^2) - 7/10).^2 ...
	+ P(:,3).^2 - 9/100;
range = [-1,1;-1,1;-0.31,0.31]; %bounds of Omega as in (42)
Nx = 11; tol = 1e-10; xcut = -0.5;
[P,h,BulkElements,SurfElements, ElementsPlot] =  ...
	generate_mesh3d(level_fun,range,Nx,tol,xcut);

% STEP 2: Matrix assembly
[K,C,M,KS,CS,MS,R]=assembly3d(P,BulkElements,SurfElements);

% STEP 3: Set model and solve PDE
a=0.1; b=0.9; alpha1=5/12; alpha2=5; beta1=5/12; beta2=0;
kappa1 = 0; kappa2 = 5; dOmega = [1;10]; dGamma = [1;10];
gammaOmega = 300; gammaGamma = 300;
f1 = @(u) a - u(:,1) + u(:,1).^2.*u(:,2);
f2 = @(u) b - u(:,1).^2.*u(:,2);
h1 = @(u,v) alpha1*v(:,1) - beta1*u(:,1) - kappa1*u(:,2);
h2 = @(u,v) alpha2*v(:,1) - beta2*u(:,1) - kappa2*u(:,2);
f = [{@(u,P,t) gammaOmega*f1(u)}; ...
     {@(u,P,t) gammaOmega*f2(u)}];
g = [{@(u,v,P,t) gammaGamma*(f1(u) - h1(u,v))}; ...
     {@(u,v,P,t) gammaGamma*(f2(u) - h2(u,v))}];
h = [{@(u,v,P,t) h1(u,v)}; {@(u,v,P,t) h2(u,v)}];
rng(0); % initialize random number generator
u0 = [a + b + 1e-3*(2*rand(size(P,1),1)-1), ...
      b/(a+b)^2 + 1e-3*(2*rand(size(P,1),1)-1)];
v0 = [a + b + 1e-3*(2*rand(size(R,2),1)-1), ...
      b/(a+b)^2 + 1e-3*(2*rand(size(R,2),1)-1)];
T = 5; tau = 1e-5;
[u,v] = solver_parabolic_bulk_surf(dOmega, dGamma, f, g, ...
	h, P, M, MS, K, KS, R, T, tau, u0, v0);
	
% STEP 4: Plot all components of numerical solution
figure, set(gcf, 'color', 'white')
sgtitle('Bulk-surface reaction-diffusion system on the torus')
subplot(2,2,1) % plot bulk Component u
plot_bulk3d(ElementsPlot, u(:,1), '$u_1$');
subplot(2,2,2) % plot bulk Component v
plot_bulk3d(ElementsPlot, u(:,2), '$u_2$')
subplot(2,2,3) % plot surface component r
plot_surf3d(P,R,SurfElements,v(:,1),'$v_1$')
xlim([xcut, 1]) % to cut the surface
subplot(2,2,4) % plot surface Component s
plot_surf3d(P,R,SurfElements,v(:,2),'$v_2$')
xlim([xcut, 1]) % to cut the surface