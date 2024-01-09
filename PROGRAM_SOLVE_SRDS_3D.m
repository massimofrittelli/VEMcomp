% DESCRIPTION - Solves SRDS on the 3D 
% unit sphere using SFEM on a triangulated mesh and plots the solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% v1_t -         LaplaceBeltrami(u1) = gamma_Gamma*g1(v1,v2),    x in \Gamma
% v2_t - d_Gamma*LaplaceBeltrami(v2) = gamma_Gamma*g2(v1,v2),    x in \Gamma
%
% \Omega = unit sphere
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


close all

% Set model parameters
a = 0.1;
b = 0.9;
d_Gamma = 10;
gamma_Gamma = 300;

% Set time discretisation
T = 5;
tau = 1e-5;

% Generate mesh
% level_fun = @(P) P(:,1).^2/4+P(:,2).^2/1+P(:,3).^2/1.5^2-1;
% [P,SurfElements] = distmeshsurface(level_fun,@huniform,0.1,...
% 	[-2.1,-1.1,-1.6; 2.1,1.1,1.6]);
xcut = -0.5;

% Matrix assembly
% [~, ~, ~, KS, MS, ~, R] = assembly_3d(P, [], SurfElements);

g1 = @(v1,v2) a - v1 + v1.^2.*v2;
g2 = @(v1,v2) b - v1.^2.*v2;

tic
NT = ceil(T/tau);
rng(0);
v1 =  a+b + 1e-3*(2*rand(length(P),1)-1);
v2 =  b/(a+b)^2 + 1e-3*(2*rand(length(P),1)-1);

MIT1 = MS+tau*KS;
MIT2 = MS+tau*d_Gamma*KS;
perm1 = symamd(MIT1);
perm2 = symamd(MIT2);
[L1,U1] = lu(MIT1(perm1,perm1),'vector');
[L2,U2] = lu(MIT2(perm2,perm2),'vector');

for i=0:NT-1

   F1 = MS*(v1 + tau*gamma_Gamma*g1(v1,v2));
   F2 = MS*(v2 + tau*gamma_Gamma*g2(v1,v2));

   v1(perm1) = U1\(L1\F1(perm1));
   v2(perm2) = U2\(L2\F2(perm2));

   if sum(isnan(v1))
        error('NAN!')
   end

end
toc

% % Plotting Numerical Solution
figure
set(gcf, 'color', 'white')

xline = ones(200,1)*xcut;
yline = linspace(-sqrt(1-xcut^2/4), sqrt(1-xcut^2/4), 200);
zline = real(1.5*sqrt(1-xcut^2/4-yline.^2));

% % Plotting Numerical Solution - Surface Component v1
subplot(1,2,1)
% Below code makes more transparent edges on the surf
trisurf(SurfElements, P(:,1), P(:,2), P(:,3), R*v1, 'FaceColor', 'interp', 'EdgeAlpha',.3)
title('$v_1$', 'Interpreter','latex')
xlabel('x')
ylabel('y')
zlabel('z', 'Rotation',0)
set(gca, 'FontSize', 18)
% plot_surf_3d(P,speye(length(P)),SurfElements,v1,'$v_1$')
% hold on
% plot3(xline,yline,zline,'k','LineWidth',1)
% plot3(xline,yline,-zline,'k','LineWidth',1)
axis equal
%xlim([xcut, max(P(:,1))])

% % Plotting Numerical Solution - Surface Component v2
subplot(1,2,2)
% Below code makes more transparent edges on the surf
trisurf(SurfElements, P(:,1), P(:,2), P(:,3), R*v2, 'FaceColor', 'interp', 'EdgeAlpha',.3)
title('$v_2$', 'Interpreter','latex')
xlabel('x')
ylabel('y')
zlabel('z', 'Rotation',0)
set(gca, 'FontSize', 18)
% plot_surf_3d(P,speye(length(P)),SurfElements,v2,'$v_2$')
% hold on
% plot3(xline,yline,zline,'k','LineWidth',1)
% plot3(xline,yline,-zline,'k','LineWidth',1)
axis equal
%xlim([xcut, max(P(:,1))])





% PLOT FOR TETRAHEDRAL MESHES

% EScut = zeros(0,3);
% for i=1:length(EB)
%     if max(PB(EB(i,:),1) >= -0.4) && max(PB(EB(i,:),1) < -0.4)
%        EScut = [EScut; EB(i,[1 2 3]); EB(i,[1 2 4]); EB(i,[1 3 4]); EB(i,[2 3 4])]; %#ok
%     end
% end
% picked = false(length(EScut),1);
% for i=1:length(EScut)
%     picked(i) = max(vecnorm(PB(EScut(i,:),:),2,2) < 0.99);
% end
% EScut = EScut(picked,:);
% 
% picked = false(length(ES),1);
% for i=1:length(ES)
%     picked(i) = max(PB(ES(i,:),1) >= -0.4);
% end
% EScut_exterior = ES(picked,:);
% 
% % Plotting Numerical Solution
% figure
% 
% % Bulk Component u
% subplot(2,2,1)
% set(gcf, 'Color','white','renderer','zbuffer')
% trisurf(EScut, PB(:,1), PB(:,2), PB(:,3), u, 'FaceColor', 'interp', 'EdgeColor', 'interp');
% hold on
% rho = 1.01;
% trisurf(EScut_exterior, rho*PB(:,1), rho*PB(:,2), rho*PB(:,3), u, 'FaceColor', 'interp', 'EdgeColor', 'none');
% set(gca,'FontSize',18)
% xlabel('x')
% ylabel('y')
% zlabel('z','rot',0)
% title('$U$', 'interpreter', 'latex')
% colorbar
% colormap jet
% axis equal
% lightangle(gca,-40,0)
% lighting gouraud
% lighting gouraud
% 
% % Bulk Component v
% subplot(2,2,2)
% trisurf(EScut, PB(:,1), PB(:,2), PB(:,3), v, 'FaceColor', 'interp', 'EdgeColor', 'interp');
% hold on
% rho = 1.01;
% trisurf(EScut_exterior, rho*PB(:,1), rho*PB(:,2), rho*PB(:,3), v, 'FaceColor', 'interp', 'EdgeColor', 'none');
% set(gca,'FontSize',18)
% xlabel('x')
% ylabel('y')
% zlabel('z','rot',0)
% title('$V$', 'interpreter', 'latex')
% colorbar
% colormap jet
% axis equal
% lightangle(gca,-40,0)
% lighting gouraud
% lighting gouraud
% 
% % Surface Component r
% subplot(2,2,3)
% trisurf(ES, PB(:,1), PB(:,2), PB(:,3), R*r, 'EdgeColor', 'none', 'FaceColor', 'interp')
% view(3)
% set(gca,'FontSize',18)
% xlabel('x')
% ylabel('y')
% zlabel('z','rot',0)
% axis equal
% xlim([-0.5,1])
% title('$R$', 'interpreter', 'latex')
% colorbar
% colormap jet
% caxis([min(r), max(r)])
% hold on
% Ccirc = [-0.5, 0, 0];   % Center of circle 
% Rcirc = sqrt(3)/2;      % Radius of circle 
% teta = 0:0.01:2*pi;
% xcirc = Ccirc(1) + zeros(size(teta)) ;
% ycirc = Ccirc(2) + Rcirc*cos(teta);
% zcirc = Ccirc(3) + Rcirc*sin(teta) ;
% plot3(xcirc,ycirc,zcirc,'k','LineWidth',1.5)
% lightangle(gca,-40,0)
% lighting gouraud
% lighting gouraud
% 
% 
% % Plotting Numerical Solution - Surface Component s
% subplot(2,2,4)
% trisurf(ES, PB(:,1), PB(:,2), PB(:,3), R*s, 'EdgeColor', 'none', 'FaceColor', 'interp')
% view(3)
% set(gca,'FontSize',18)
% xlabel('x')
% ylabel('y')
% zlabel('z','rot',0)
% axis equal
% xlim([-0.5,1])
% title('$S$', 'interpreter', 'latex')
% colorbar
% colormap jet
% caxis([min(s), max(s)])
% hold on
% Ccirc = [-0.5, 0, 0];   % Center of circle 
% Rcirc = sqrt(3)/2;      % Radius of circle 
% teta = 0:0.01:2*pi;
% xcirc = Ccirc(1) + zeros(size(teta)) ;
% ycirc = Ccirc(2) + Rcirc*cos(teta);
% zcirc = Ccirc(3) + Rcirc*sin(teta) ;
% plot3(xcirc,ycirc,zcirc,'k','LineWidth',1.5)
% lightangle(gca,-40,0)
% lighting gouraud
% lighting gouraud