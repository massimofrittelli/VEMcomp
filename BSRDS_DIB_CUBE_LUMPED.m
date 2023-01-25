% DESCRIPTION - Solves the BS DIB model on the 3D cube with VEM on a cubic
% mesh.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bulk-surface DIB model in 3D.
%
% Notations:
% Omega = [0,L]^3               (cube)
% F_0 = [0,L]^2 x {0}           (bottom face)
% F_1 = [0,L]^2 x {L}           (top face)
% Gamma = \partial \Omega \setminus (F_0 \cup F_1)          (lateral faces)
%
%
% Bulk system in Omega:
% b_t     =         Laplace(b)              - k_b*(b-b_bulk),       in \Omega
% q_t     =  d_bulk*Laplace(q)              - k_q*(q-q_bulk),       in \Omega
%
% Surface system on F_0:
% eta_t   =         Laplace_Beltrami(eta)   + f_3(b,eta,theta),     on F_0
% theta_t = d_gamma*Laplace_Beltrami(theta) + f_4(q,eta,theta),     on F_0
% 
% Coupling boundary conditions on F_0:
% flux(b) = -f_3(b,eta,theta)*psi_eta,                              on F_0
% flux(q) = -f_4(q,eta,theta)*psi_theta,                            on F_0
%
% Dirichlet boundary conditions on F_1:
% b = b_bulk                                                        on F_1
% q = q_bulk                                                        on F_1
%
% Neumann boundary conditions on \Gamma
% flux(b) = 0,                                                      on \Gamma
% flux(q) = 0,                                                      on \Gamma
%
% Kinetics definitions:
% f_3 := rho * [b*A_1*(1-theta)*eta -A_2*eta^3 -B(theta-alpha)]
% f_4 := rho * [C*q*(1+k_2*eta)*(1-theta)*(1-gamma*(1-theta)) - D*(1+k_3*eta)*theta*(1+gamma*theta)] 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       BEGIN INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% LOAD MESH
load('cube33.mat')
L = 10;         % cube edge length
P = P*L;        % rescaled  nodes
% Non-Lumped version
% M = M*L^3;      % rescaled bulk mass matrix
% Mbot = Mbot*L^2;    % rescaled surface mass matrix
% Lumped version
M = spdiags(sum(M*L^3,2),0,length(M),length(M));      % rescaled bulk mass matrix
Mbot = spdiags(sum(Mbot*L^2,2),0,length(Mbot),length(Mbot));    % rescaled surface mass matrix
K = K*L;        % rescaled bulk stiffness matrix
% no rescaling needed for surface stiffness matrix

% SET FINAL TIME AND TIME STEP
T = 100; % 100
tau = 5e-3;

% SET DIFFUSION COEFFICIENTS
d_bulk = 1;
d_gamma = 20;

% SET BULK REACTION KINETICS
b_bulk = 1;
q_bulk = 1;
k_b = 1;
k_q = 1;

% SET DIB REACTION KINETICS
rho = 1.5; %1
A_1 = 10;
A_2 = 1;
alpha = 0.5;
k2 = 2.5;
k3 = 1.5;
B = 30;
C = 3;
gamma = .2;
D = C*(1-alpha)*(1-gamma+gamma*alpha)/(alpha*(1+gamma*alpha));

% FIGURA 4b EJAM
% FIGURA 2 AMOS

% IL SEGUENTE E' UN PARAMETRO DI PERTURBAZIONE DEL PROBLEMA B-S RISPETTO AL
% PROBLEMA 1D
pert_param = 1.1e-1;

% SET BULK-SURFACE COUPLING PARAMETERS
psi_eta = pert_param;
psi_theta = pert_param;

% SET INITIAL CONDITIONS
rng(0)
ampl = 1e-2;
b0 = @(x) 0*pert_param*b_bulk/L*x(:,2) + (1-0*pert_param)*b_bulk;
q0 = @(x) 0*pert_param*b_bulk/L*x(:,2) + (1-0*pert_param)*b_bulk;
eta0 = @(x) ampl*(1+rand(size(x(:,1))))/2;
theta0 = @(x) ampl*rand(size(x(:,1))) + alpha;

% SET PLOT FREQUENCY (ONE FRAME PER 'freq' ITERATIONS)
plot_freq = 100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       END INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% To enforce Dirichlet boundary conditions, we define auxiliary variables
% b_dir := b - b_bulk and q_dir := q - q_bulk and we solve the following 
% auxiliary problem 
%
% Bulk system in Omega:
% b_dir_t     =         Laplace(b_dir)      - k_b*(b_dir),          in \Omega
% q_dir_t     =  d_bulk*Laplace(q_dir)      - k_q*(q_dir),          in \Omega
%
% Surface system on F_0:
% eta_t   =         Laplace_Beltrami(eta)   + f_3(b_dir+b_bulk,eta,theta),     on F_0
% theta_t = d_gamma*Laplace_Beltrami(theta) + f_4(q_dir+q_bulk,eta,theta),     on F_0
% 
% Coupling boundary conditions on F_0:
% flux(b_dir) = -f_3(b_dir+b_bulk,eta,theta)*psi_eta,               on F_0
% flux(q_dir) = -f_4(q_dir+q_bulk,eta,theta)*psi_theta,             on F_0
%
% Dirichlet boundary conditions on F_1:
% b_dir = 0                                                         on F_1
% q_dir = 0                                                         on F_1
%
% Neumann boundary conditions on \Gamma
% flux(b_dir) = 0,                                                  on \Gamma
% flux(q_dir) = 0,                                                  on \Gamma
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% BULK KINETICS
f_1 = @(bdir) -k_b*bdir;
f_2 = @(qdir) -k_q*qdir;

% SURFACE KINETICS (MODIFIED DIB KINETICS)
f_3 = @(bdir, eta, theta) rho*(A_1*(Abot'*Atop*bdir+b_bulk).*(1-theta).*eta - A_2*eta.^3 - B*(theta-alpha));
f_4 = @(qdir, eta, theta) rho*(C*(Abot'*Atop*qdir+q_bulk).*(1+k2*eta).*(1-theta).*(1-gamma*(1-theta)) - D*(1+k3*eta).*theta.*(1+gamma*theta));

% PURE DIB KINETICS
f_eta = @(eta, theta) rho*(A_1*(1-theta).*eta - A_2*eta.^3 - B*(theta-alpha));
f_theta = @(eta, theta) rho*(C*(1+k2*eta).*(1-theta).*(1-gamma*(1-theta)) - D*(1+k3*eta).*theta.*(1+gamma*theta));

% FLUX TERMS
flux1 = @(bdir, eta, theta) -f_3(bdir, eta, theta)*psi_eta;
flux2 = @(qdir, eta, theta) -f_4(qdir, eta, theta)*psi_theta;

% AUXILIARY MATRICES FOR TIME STEPPING
MIT1  = (Atop'*(M+tau*K)*Atop);
MIT1d = (Atop'*(M+d_bulk*tau*K)*Atop);
MIT3  = (Mbot+tau*Kbot);
MIT3d = (Mbot+d_gamma*tau*Kbot);
perm1 = symamd(MIT1);
perm1d = symamd(MIT1d);
perm3 = symamd(MIT3);
perm3d = symamd(MIT3d);
[L1, U1]  = lu(MIT1(perm1,perm1),'vector');
[L1d,U1d] = lu(MIT1d(perm1d,perm1d),'vector');
[L3, U3]  = lu(MIT3(perm3,perm3),'vector');
[L3d,U3d] = lu(MIT3d(perm3d,perm3d),'vector');

M1 = (Atop'*M*Atop);
M2 = (Atop'*Abot*Mbot);
M3 = Mbot;

% INITIALIZING NUMERICAL SOLUTION
N = size(Atop,2);
n = size(Abot,2);
bdir = b0(Atop'*P) - b_bulk;
qdir = q0(Atop'*P) - q_bulk;
eta = eta0(Abot'*P);
theta = theta0(Abot'*P);
eta_pure = eta;
theta_pure = theta;

% PLOTTING INITIAL CONDITIONS
cut = P(:,1) <= L/2;
Pcut = P(cut,:);
TRIcut = convhull(Pcut);
TRI = delaunay(Abot'*P(:,1), Abot'*P(:,2));

figure(1)
set(gcf, 'Renderer', 'zbuffer', 'color','white','Position', [100 100 500 400])
b = Atop*bdir + b_bulk;
trisurf(TRIcut, Pcut(:,1), Pcut(:,2), Pcut(:,3), b(cut),'EdgeColor','none','FaceColor', 'interp');
hold on
trisurf(TRI, Abot'*P(:,1), Abot'*P(:,2), Abot'*P(:,3), eta,'EdgeColor','none','FaceColor', 'interp');
colorbar('FontSize', 18)
colormap jet
xlim([0,L])
ylim([0,L])
zlim([0,L])
set(gca, 'fontsize',18)
title('$\eta, b$', 'interpreter', 'latex')
caxis([min([min(b),min(eta),min(eta_pure)]), max([max(b),max(eta),max(eta_pure)])])

figure(3)
set(gcf, 'Renderer', 'zbuffer', 'color','white','Position', [100 100 500 400])
trisurf(TRI, Abot'*P(:,1), Abot'*P(:,2), Abot'*P(:,3), eta_pure,'EdgeColor','none','FaceColor', 'interp');
xlim([0,L])
ylim([0,L])
zlim([0,L])
set(gca, 'fontsize',18)
caxis([min([min(bdir+b_bulk),min(eta),min(eta_pure)]), max([max(bdir+b_bulk),max(eta),max(eta_pure)])])
colorbar('FontSize', 18)
title('$\eta$', 'interpreter', 'latex')
colormap jet
view(2)
    
% for i=1:10
%     MOV(i) = getframe(gcf);
% end

tic
%COMPUTING NUMERICAL SOLUTION
for i=1:ceil(T/tau)
    Fb = M1*(bdir+tau*f_1(bdir)) + M2*tau*flux1(bdir,eta,theta);
    Fq = M1*(qdir+tau*f_2(qdir)) + M2*tau*flux2(qdir,eta,theta);
    Feta = M3*(eta+ tau*f_3(bdir,eta,theta));
    Ftheta = M3*(theta+ tau*f_4(qdir,eta,theta));
    Feta_pure = M3*(eta_pure + tau*f_eta(eta_pure, theta_pure));
    Ftheta_pure = M3*(theta_pure + tau*f_theta(eta_pure, theta_pure));
    bdir(perm1) = U1\(L1\Fb(perm1));
    qdir(perm1d) = U1d\(L1d\Fq(perm1d));
    eta(perm3) = U3\(L3\Feta(perm3));
    theta(perm3d) = U3d\(L3d\Ftheta(perm3d));
    eta_pure(perm3) = U3\(L3\Feta_pure(perm3));
    theta_pure(perm3d) = U3d\(L3d\Ftheta_pure(perm3d));

    
    if rem(i,plot_freq) == 0
        
        fprintf('t=%d\n', tau*i)
        
        % PLOTTING b,eta OF COUPLED BSRDS MODEL
        figure(1)
        cla
        b = Atop*bdir + b_bulk;
        trisurf(TRIcut, Pcut(:,1), Pcut(:,2), Pcut(:,3), b(cut),'EdgeColor','none','FaceColor', 'interp');
        hold on
        trisurf(TRI, Abot'*P(:,1), Abot'*P(:,2), Abot'*P(:,3), eta,'EdgeColor','none','FaceColor', 'interp');
        xlim([0,L])
        ylim([0,L])
        zlim([0,L])
        set(gca, 'fontsize',18)
        colorbar('FontSize', 18)
        title('$\eta, b$', 'interpreter', 'latex')
        caxis([min([min(b),min(eta),min(eta_pure)]), max([max(b),max(eta),max(eta_pure)])])
        
        % PLOTTING eta OF UNCOUPLED SRDS MODEL
        figure(3)
        cla
        trisurf(TRI, Abot'*P(:,1), Abot'*P(:,2), Abot'*P(:,3), eta_pure,'EdgeColor','none','FaceColor', 'interp');
        xlim([0,L])
        ylim([0,L])
        set(gca, 'fontsize',18)
        colorbar('FontSize', 18)
        title('$\eta$', 'interpreter', 'latex')
        caxis([min([min(bdir+b_bulk),min(eta),min(eta_pure)]), max([max(bdir+b_bulk),max(eta),max(eta_pure)])])
        view(2)
    end
end
toc

% PLOTTING b,eta OF COUPLED BSRDS MODEL AT FINAL TIME WITN NON_SCALED COLORMAPS
figure(4)
set(gcf, 'Renderer', 'zbuffer', 'color','white','Position', [100 100 500 400])
b = Atop*bdir + b_bulk;
b = b(cut)*(max([max(eta),max(eta_pure)])-min([min(eta),min(eta_pure)]))/(max(b(cut))-min(b(cut)));
b = b-min(b) + min([min(eta),min(eta_pure)]);
trisurf(TRIcut, Pcut(:,1), Pcut(:,2), Pcut(:,3), b,'EdgeColor','none','FaceColor', 'interp');
hold on
trisurf(TRI, Abot'*P(:,1), Abot'*P(:,2), Abot'*P(:,3), eta,'EdgeColor','none','FaceColor', 'interp');
xlim([0,L])
ylim([0,L])
zlim([0,L])
set(gca, 'fontsize',18)
title('$\eta, b$', 'interpreter', 'latex')
caxis([min([min(eta),min(eta_pure)]), max([max(eta),max(eta_pure)])])
colormap jet  
colorbar off

% movie2avi(MOV, 'BS_DIB_3D_cube.avi')