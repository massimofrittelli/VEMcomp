%% VEMCOMP:  Example 6.2.2 surface-only reaction-diffusion system (eq. 33)
%% NOTICE: Add DistMesh to the current folder before running this script
% Reaction-diffusion system on an ellipsoidal surface, where the mesh is
% a triangulated surface generated using DistMesh
clear all
disp('Reaction-diffusion system on ellipsoidal surface')

% STEP 1: Generate mesh using DistMesh
h = 0.1; % Upper bound of the meshsize
level_fun = @(P) P(:,1).^2/4+P(:,2).^2/1+P(:,3).^2/1.5^2-1;
[P,SurfElements] = distmeshsurface(level_fun,@huniform,h,...
	[-2.1,-1.1,-1.6; 2.1,1.1,1.6]);
	
% STEP 2: Matrix assembly
[~, ~, ~, KS, MS,~,R] = assembly3d(P, [], SurfElements);
% There are no bulk elements and no bulk matrices

% STEP 3: Set model and solve PDE
a = 0.1; b = 0.9; dGamma = [1;10]; gammaGamma = 300;
g = [{@(u,P,t) gammaGamma*(a-u(:,1)+u(:,1).^2.*u(:,2))}; ...
     {@(u,P,t) gammaGamma*(b-u(:,1).^2.*u(:,2))}];
rng(0); % initialize random number generator
v0 = [a + b + 1e-3*(2*rand(size(P,1),1)-1), ...
      b/(a+b)^2 + 1e-3*(2*rand(size(P,1),1)-1)];
T = 10; tau = 1e-5;
[v, t, vprime_norm, v_average] = solver_parabolic_surf(dGamma,g,P,MS,KS,R,T,tau,v0);

% STEP 4: Post-processing
figure, set(gcf, 'color', 'white')
sgtitle('Reaction-diffusion system on ellipsoidal surface')
subplot(2,2,1) % plot component v1
plot_surf3d(P,R,SurfElements,v(:,1),'$v_1$')
subplot(2,2,2) % plot component v2
plot_surf3d(P,R,SurfElements,v(:,2),'$v_2$')
subplot(2,2,3), semilogy(t, vprime_norm, 'LineWidth',2)
legend('$\|\dot{v}_1\|_{L^2(\Gamma_h)}$', 'interpreter', 'latex') 
xlabel('t'), set(gca, 'FontSize', 14)
subplot(2,2,4), plot(t, v_average, 'LineWidth',2)
legend('$<v_1>$', 'interpreter', 'latex') 
xlabel('t'), set(gca, 'FontSize', 14)

