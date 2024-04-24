% DESCRIPTION - Solves Anotida's 4-species Schnakenberg BSRDS on the 3D 
% unit sphere using VEM on a polyhedral mesh and plots the solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% u_t -         Laplace(u)         = gamma_Omega f(u,v)                    x in \Omega
% v_t - d_Omega Laplace(v)         = gamma_Omega g(u,v)                    x in \Omega
% r_t -         LaplaceBeltrami(r) = gamma_Gamma(f(r,s) - h1(u,v,r,s)),    x in \Gamma
% s_t - d_Gamma LaplaceBeltrami(s) = gamma_Gamma(g(r,s) - h2(u,v,r,s)),    x in \Gamma
%         \nabla u \cdot \nu = gamma_\Gamma h1(u,v,r,s),                   x in \Gamma
% d_Omega \nabla v \cdot \nu = gamma_\Gamma h2(u,v,r,s),                   x in \Gamma
%
% \Omega = unit sphere
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


close all
clearvars

% Set model parameters
a = 0.1;
b = 0.9;
alpha1 = 5/12;
alpha2 = 5;
beta1 = 5/12;
beta2 = 0;
k1 = 0;
k2 = 5;

d_Omega = 10;
d_Gamma = 10;
gamma_Omega = 300;
gamma_Gamma = 300;

% Set time discretisation
T = 1e-5;
tau = 1e-5;

% Set space discretisation
% Generating mesh
fun = @(P) P(:,1).^2 + P(:,2).^2 + P(:,3).^2 - 1;
range = [-1,1; -1,1; -1,1];
tol = 1e-3;
Nx = 10;
xcut = -0.5;
[P, h, BulkElements, SurfElements, ElementsPlot] = generate_mesh3d(fun, range, Nx, tol, xcut);

% Assembling matrices
[K, M, C, KS, MS, CS, R] = assembly3d(P, BulkElements, SurfElements);

N = length(P);
NGamma = length(MS); % Amount of boundary nodes

f = @(u,v) a - u + u.^2.*v;
g = @(u,v) b - u.^2.*v;
h1 = @(u,v,r) alpha1*r - beta1*u - k1*v;
h2 = @(u,v,s) alpha2*s - beta2*u - k2*v;

tic
NT = ceil(T/tau);
rng(0);
u =  a+b + 1e-3*(2*rand(length(P),1)-1);
v =  b/(a+b)^2 + 1e-3*(2*rand(length(P),1)-1);
r =  a+b + 1e-3*(2*rand(size(R,2),1)-1);
s =  b/(a+b)^2 + 1e-3*(2*rand(size(R,2),1)-1);

MIT1 = M+tau*K;
MIT2 = M+tau*d_Omega*K;
MIT3 = MS+tau*KS;
MIT4 = MS+tau*d_Gamma*KS;
perm1 = symamd(MIT1);
perm2 = symamd(MIT2);
perm3 = symamd(MIT3);
perm4 = symamd(MIT4);
[L1,U1] = lu(MIT1(perm1,perm1),'vector');
[L2,U2] = lu(MIT2(perm2,perm2),'vector');
[L3,U3] = lu(MIT3(perm3,perm3),'vector');
[L4,U4] = lu(MIT4(perm4,perm4),'vector');

for i=0:NT-1
   
%    unext = (M+tau*K)          \(M* (u + tau*gamma_Omega*f(u,v)) + tau*gamma_Gamma*R*MS*h1(R'*u,R'*v,r));
%    vnext = (M+tau*d_Omega*K)  \(M* (v + tau*gamma_Omega*g(u,v)) + tau*gamma_Gamma*R*MS*h2(R'*u,R'*v,s));
%    rnext = (MS+tau*KS)        \(MS*(r + tau*gamma_Gamma*(f(r,s) - h1(R'*u,R'*v,r))));
%    s     = (MS+tau*d_Gamma*KS)\(MS*(s + tau*gamma_Gamma*(g(r,s) - h2(R'*u,R'*v,s))));

%    u = unext;
%    v = vnext;
%    r = rnext;

   F1 = M* (u + tau*gamma_Omega*f(u,v)) + tau*gamma_Gamma*R*MS*h1(R'*u,R'*v,r);
   F2 = M* (v + tau*gamma_Omega*g(u,v)) + tau*gamma_Gamma*R*MS*h2(R'*u,R'*v,s);
   F3 = MS*(r + tau*gamma_Gamma*(f(r,s) - h1(R'*u,R'*v,r)));
   F4 = MS*(s + tau*gamma_Gamma*(g(r,s) - h2(R'*u,R'*v,s)));

   u(perm1) = U1\(L1\F1(perm1));
   v(perm2) = U2\(L2\F2(perm2));
   r(perm3) = U3\(L3\F3(perm3));
   s(perm4) = U4\(L4\F4(perm4));

   if sum(isnan(r))
        error('NAN!')
   end

end
toc

% % Plotting Numerical Solution
figure
set(gcf, 'color', 'white')

% % Plotting Numerical Solution - Bulk Component u
subplot(2,2,1)
plot_bulk3d(ElementsPlot, u, '$u$')
% Below code makes thicker lines on the cut
% for i=1:size(ElementsPlot)
%     hold on
%     if not(ElementsPlot(i).is_boundary)
%         plot(ElementsPlot(i), u(ElementsPlot(i).Pind));
%     else
%         plot(ElementsPlot(i), u(ElementsPlot(i).Pind),'k', 0.3);
%     end
% end


% % Plotting Numerical Solution - Bulk Component v
subplot(2,2,2)
plot_bulk3d(ElementsPlot, v, '$v$')
% Below code makes thicker lines on the cut
% for i=1:size(ElementsPlot)
%     hold on
%     if not(ElementsPlot(i).is_boundary)
%         plot(ElementsPlot(i), v(ElementsPlot(i).Pind));
%     else
%         plot(ElementsPlot(i), v(ElementsPlot(i).Pind),'k', 0.3);
%     end
% end


xline = ones(200,1)*xcut;
yline = linspace(-sqrt(3)/2, sqrt(3)/2, 200);
zline = sqrt(9/100-(sqrt(xcut^2+yline.^2)-7/10).^2);


% % Plotting Numerical Solution - Surface Component r
subplot(2,2,3)
plot_surf3d(P,R,SurfElements,r,'$r$')
hold on
plot3(xline,yline,zline,'k','LineWidth',1)
plot3(xline,yline,-zline,'k','LineWidth',1)
xlim([xcut, max(P(:,1))])
% Below code makes more transparent esdes on the surf
% trisurf(SurfaceElements, P(:,1), P(:,2), P(:,3), R*r, 'FaceColor', 'interp', 'EdgeAlpha',.3)


% % Plotting Numerical Solution - Surface Component s
subplot(2,2,4)
plot_surf3d(P,R,SurfElements,s,'$s$')
hold on
plot3(xline,yline,zline,'k','LineWidth',1)
plot3(xline,yline,-zline,'k','LineWidth',1)
xlim([xcut, max(P(:,1))])
% Below code makes more transparent esdes on the surf
% trisurf(SurfaceElements, P(:,1), P(:,2), P(:,3), R*s, 'FaceColor', 'interp', 'EdgeAlpha',.3)





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