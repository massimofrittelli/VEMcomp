clear variables
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bulk-surface DIB model in 3D.
%
% Notations:
% Omega = [0,L]^3               (cube)
% F_0 = [0,L]^2 x {0}           (bottom face)
% F_1 = [0,L]^2 x {1}           (top face)
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
load('cube6.mat')
L = 10;
P = P*L;
M = M*L^3;
MS = MS*L^2;
K = K*L;
%KS = KS; % No transformation needed

% SET FINAL TIME AND TIME STEP
T = 150;
tau = 1e-2;

% SET DIFFUSION COEFFICIENTS
d_bulk = 1;
d_gamma = 20;

% SET BULK REACTION KINETICS
b_bulk = 1;
q_bulk = 1;
k_b = 1;
k_q = 1;

% SET DIB REACTION KINETICS
rho = 1;
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
M1 = (Atop'*(M+tau*K)*Atop)\(Atop'*M*Atop);
M1d = (Atop'*(M+d_bulk*tau*K)*Atop)\(Atop'*M*Atop);

M2 = (Atop'*(M+tau*K)*Atop)\(Atop'*Abot*MS);
M2d = (Atop'*(M+d_bulk*tau*K)*Atop)\(Atop'*Abot*MS);

M3 = (MS+tau*KS)\MS;
M3d = (MS+d_gamma*tau*KS)\MS;



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
figure('Renderer', 'zbuffer', 'color','white','Position', [100 100 500 400])
xlim([0,L])
ylim([0,L])
zlim([0,L])
scatter3(P(:,1), P(:,2), P(:,3), Atop*bdir+b_bulk);
colorbar('FontSize', 18)
caxis([min([min(bdir+b_bulk),min(eta),min(eta_pure)]), max([max(bdir+b_bulk),max(eta),max(eta_pure)])])

% TODO: plot correctly surf solution (mesh probably needed, or use Delaunay command)
xlim([0,L])
ylim([0,L])
zlim([0,L])
plot_surface(Abot'*P, "none", eta)
colorbar('FontSize', 18)
caxis([min([min(bdir+b_bulk),min(eta),min(eta_pure)]), max([max(bdir+b_bulk),max(eta),max(eta_pure)])])

figure('Renderer', 'zbuffer', 'color','white','Position', [100 100 500 400])
set(gca, 'fontsize',18)
xlim([0,L])
ylim([0,L])
zlim([0,L])
plot_surface(Abot'*P, "none", eta_pure)
caxis([min([min(bdir+b_bulk),min(eta),min(eta_pure)]), max([max(bdir+b_bulk),max(eta),max(eta_pure)])])
colorbar('FontSize', 18)
    
% for i=1:10
%     MOV(i) = getframe(gcf);
% end


%COMPUTING NUMERICAL SOLUTION
for i=1:ceil(T/tau)
    bnew = M1*(bdir+tau*f_1(bdir)) + M2*tau*flux1(bdir,eta,theta);
    qnew = M1d*(qdir+tau*f_2(qdir)) + M2d*tau*flux2(qdir,eta,theta);
    etanew = M3*(eta+ tau*f_3(bdir,eta,theta));
    thetanew = M3d*(theta+ tau*f_4(qdir,eta,theta));
    eta_pure_new = M3*(eta_pure + tau*f_eta(eta_pure, theta_pure));
    theta_pure_new = M3d*(theta_pure + tau*f_theta(eta_pure, theta_pure));
    bdir = bnew;
    qdir = qnew;
    eta = etanew;
    theta = thetanew;
    eta_pure = eta_pure_new;
    theta_pure = theta_pure_new;

    
    if rem(i,plot_freq) == 0
        figure(1)
        cla
        plot_bulk_polygons(P,E,'b','b',Atop*bdir+b_bulk);
        caxis([min([min(bdir+b_bulk),min(eta),min(eta_pure)]), max([max(bdir+b_bulk),max(eta),max(eta_pure)])])
        plot_surface(Abot'*P, "none", eta)
        caxis([min([min(bdir+b_bulk),min(eta),min(eta_pure)]), max([max(bdir+b_bulk),max(eta),max(eta_pure)])])
        %MOV(i/plot_freq+10) = getframe(gcf);
        
        figure(2)
        cla
        xlabel('x','FontSize',18)
        ylabel('y','FontSize',18, 'rot', 0)
        plot_surface(Abot'*P, "none", eta_pure)
        caxis([min([min(bdir+b_bulk),min(eta),min(eta_pure)]), max([max(bdir+b_bulk),max(eta),max(eta_pure)])])
    end
end

% movie2avi(MOV, 'receptor_ligand.avi')