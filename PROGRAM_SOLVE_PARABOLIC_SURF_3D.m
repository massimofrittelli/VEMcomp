% DESCRIPTION - Solves a parabolic surface toy model on the sphere, plots
% solution and computes errors.

close all
clearvars

load('mesh_sphere_Nx11.mat')
esol_v = @(P,t) P(:,1).*P(:,2).*P(:,3)*exp(t);

T = 1;
tau = 1e-4;
g = @(P,t) 13*P(:,1).*P(:,2).*P(:,3)*exp(t);
v0 = R'*(P(:,1).*P(:,2).*P(:,3));

v = solver_parabolic3d_surf_sphere(g, R, P, MS, KS, T, tau, v0);

es_v = esol_v(R'*P,T);
err_v = v - es_v;
L2err_abs = sqrt(err_v'*MS*err_v);
normsol = sqrt(es_v'*MS*es_v);
L2err_rel = L2err_abs/normsol;


% Plotting Numerical Solution
figure
set(gcf,'Color','white')
plot_surf_3d(P,R,SurfaceElements,v,'$v$')
xlim([-0.5,1])
hold on
Ccirc = [-0.5, 0, 0];   % Center of circle 
Rcirc = sqrt(3)/2;      % Radius of circle 
teta = 0:0.01:2*pi;
xcirc = Ccirc(1) + zeros(size(teta));
ycirc = Ccirc(2) + Rcirc*cos(teta);
zcirc = Ccirc(3) + Rcirc*sin(teta);
plot3(xcirc,ycirc,zcirc,'k','LineWidth',1.5)

% Plotting Numerical Error
figure
set(gcf,'Color','white')
plot_surf_3d(P,R,SurfaceElements,esol_v(R'*P,T)-v,'$v_{exact}-v$')
xlim([-0.5,1])
hold on
Ccirc = [-0.5, 0, 0];   % Center of circle 
Rcirc = sqrt(3)/2;      % Radius of circle 
teta = 0:0.01:2*pi;
xcirc = Ccirc(1) + zeros(size(teta)) ;
ycirc = Ccirc(2) + Rcirc*cos(teta);
zcirc = Ccirc(3) + Rcirc*sin(teta) ;
plot3(xcirc,ycirc,zcirc,'k','LineWidth',1.5)