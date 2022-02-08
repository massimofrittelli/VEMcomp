% Script to test a parabolic coupled  bulk-surface problem on the 3D sphere

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% u_t - Laplace(u) = xyz exp(-t),                                   x in \Omega
% v_t - LaplaceBeltrami(v) + \nabla u \cdot \nu = 16xyz exp(-t),    x in \Gamma
% \nabla u \cdot \nu = 3xyz exp(-t),                                x in \Gamma
%
% \Omega = unit sphere
%
% Exact solution:
%
% u(x,t) =   xyz exp(- t)
% v(x,t) = 2 xyz exp(- t)
%
%
% Note: if the solution is exponentially-in-time increasing instead of
% decreasing, it causes stability issues with the method. 
% The solution blows up. The stability condition probably becomes too
% restrictive.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


close all
%clearvars

% Set free model parameters

% Dependent model parameters (DO NOT TOUCH);

% Set time discretisation
 T = 2^(-2);
% tau = 2^(-4);
% tau = 2^(-6);
% tau = 2^(-8);
% tau = 2^(-10);

% Set space discretisation
% load('sphere17.mat')

% END OF INPUT

N = length(PB);
NGamma = length(MS); % Amount of boundary nodes
R = spalloc(N, NGamma, NGamma);
R(boundarynode, 1:NGamma) = speye(NGamma); % reduction matrix

esol_u = @(P,t)   P(:,1).*P(:,2).*P(:,3)*exp(-t);
esol_v = @(P,t) 2*P(:,1).*P(:,2).*P(:,3)*exp(-t);

tic
NT = ceil(T/tau);
u =  esol_u(PB,0);
v =  esol_v(R'*PB,0);

MIT1 = M+tau*K;
MIT2 = MS+tau*KS;
perm1 = symamd(MIT1);
perm2 = symamd(MIT2);
[L1,U1] = lu(MIT1(perm1,perm1),'vector');
[L2,U2] = lu(MIT2(perm2,perm2),'vector');

for i=0:NT-1
    
   F1 = (M*(u - tau*u) + tau*R*MS*(v+R'*u));
   F2 = (MS*(v + tau*8*R'*u + 7*tau*v));
   
   u(perm1) =  U1\(L1\F1(perm1));
   v(perm2) =  U2\(L2\F2(perm2));

end
toc

es_u = esol_u(PB,T);
es_v = esol_v(R'*PB,T);
err_u = u - es_u;
err_v = v - es_v;
L2err_u = sqrt(err_u'*M*err_u);
L2err_v = sqrt(err_v'*MS*err_v);
L2err_prod = sqrt(err_u'*M*err_u + err_v'*MS*err_v);

% % Plotting Numerical Solution - Bulk Component u
figure
set(gcf, 'Color','white')
set(gcf, 'Renderer','zbuffer')
hold on
for i=1:size(ElementsCut)
    plotSolution(ElementsCut(i), u);
end
rho = 1.01;
trisurf(EGammaCut, rho*PB(:,1), rho*PB(:,2), rho*PB(:,3), u, 'FaceColor', 'interp', 'EdgeColor', 'none');
view(3)
set(gca,'FontSize',18)
xlabel('x')
ylabel('y')
zlabel('z','rot',0)
title('$U$', 'interpreter', 'latex')
colorbar
colormap jet
axis equal
axis([-1.1,1.1,-1.1,1.1,-1.1,1.1])
lightangle(gca,-40,0)
 

% % Plotting Numerical Solution - Surface Component v
figure
set(gcf,'Color','white')
trisurf(EGamma, PB(:,1), PB(:,2), PB(:,3), R*v, 'EdgeColor', 'none', 'FaceColor', 'interp')
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
lightangle(gca,-40,0)
% lighting gouraud


% % Plotting Numerical Error - Bulk Component u
% figure
% indsol = P(:,1) >= -0.5;
% chull = convhull(P(indsol,1), P(indsol,2), P(indsol,3));
% set(gcf, 'Color','white')
% set(gcf, 'Renderer','painters')
% %set(gcf, 'Renderer', 'zbuffer')
% trisurf(chull, P(indsol,1), P(indsol,2), P(indsol,3), es_u(indsol,1)-u(indsol,1),'EdgeColor', 'none', 'FaceColor', 'interp')
% view(3)
% set(gca,'FontSize',18)
% xlabel('x')
% ylabel('y')
% zlabel('z','rot',0)
% axis equal
% caxis([min(es_u-u), max(es_u-u)])
% title('$u^{-\ell}-U$', 'interpreter', 'latex')
% colorbar
% hold on
% Ccirc = [-0.5, 0, 0];   % Center of circle 
% Rcirc = sqrt(3)/2;      % Radius of circle 
% teta = 0:0.01:2*pi;
% xcirc = Ccirc(1) + zeros(size(teta)) ;
% ycirc = Ccirc(2) + Rcirc*cos(teta);
% zcirc = Ccirc(3) + Rcirc*sin(teta) ;
% plot3(xcirc,ycirc,zcirc,'k','LineWidth',1.5)
% % lightangle(gca,-65,30)
% % lighting gouraud
% 
% % Plotting Numerical Error - Surface Component v
% figure
% set(gcf,'Color','white')
% trisurf(EGamma, P(:,1), P(:,2), P(:,3), R*(es_v-v), 'EdgeColor', 'none', 'FaceColor', 'interp')
% view(3)
% set(gca,'FontSize',18)
% xlabel('x')
% ylabel('y')
% zlabel('z','rot',0)
% axis equal
% caxis([min(es_v-v), max(es_v-v)])
% xlim([-0.5,1])
% title('$v^{-\ell}-V$','interpreter', 'latex')
% colorbar
% hold on
% Ccirc = [-0.5, 0, 0];   % Center of circle 
% Rcirc = sqrt(3)/2;      % Radius of circle 
% teta = 0:0.01:2*pi;
% xcirc = Ccirc(1) + zeros(size(teta)) ;
% ycirc = Ccirc(2) + Rcirc*cos(teta);
% zcirc = Ccirc(3) + Rcirc*sin(teta) ;
% plot3(xcirc,ycirc,zcirc,'k','LineWidth',1.5)
% % lightangle(gca,-65,30)
% % lighting gouraud