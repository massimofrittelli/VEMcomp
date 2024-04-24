%% VEMCOMP:  RDS solver in the bulk
% Neumann boundary conditions
clearvars
disp('RDS solver, 2D domain, Neumann BCs')

%% STEP 1: Generate mesh
level_fun = @(P) P(:,1).^2 + P(:,2).^2 -1;
range = [-1,1; -1,1]; Nx = 40; tol = 1e-6;
[P, h, BulkElements, SurfElements] = generate_mesh2d(level_fun, range, Nx, tol);

%% STEP 2: Matrix assembly
[K,C,M,KS,MS,R] = assembly2d(P, BulkElements, SurfElements);

%% STEP3: Solve PDE
D = [1;10]; % diffusion coefficients
a = 0.1; b = 0.9;
rho = 300; % Effective domain area
T = 5;
tau = 1e-4;
rng(0); % Initialize seed for random number generator
u0 = [a + b + 1e-4*(2*rand(size(P,1),1)-1), ...
      b/(a+b)^2 + 1e-4*(2*rand(size(P,1),1)-1)];
f = [{@(u,P,t) rho*(a-u(:,1)+u(:,1).^2.*u(:,2))}; ...
     {@(u,P,t) rho*(b-u(:,1).^2.*u(:,2))}];
bcond = 'neu';
u = solver_parabolic_bulk(D, f , P, M, K, R, bcond ,T, tau , u0) ;

%% STEP 4: Post-processing
figure, set(gcf,'color','white')
plot_bulk2d(BulkElements, u(:,1), '$u$');
figure, set(gcf,'color','white')
plot_bulk2d(BulkElements, u(:,2), '$v$');