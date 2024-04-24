% DESCRIPTION - solves parabolic B-S problem on the 3D sphere, plots
% solution and computes error.

%
% u_t - Laplace(u) = xyz exp(-t),                                   x in \Omega
% v_t - LaplaceBeltrami(v) + \nabla u \cdot \nu = 16xyz exp(-t),    x in \Gamma
% \nabla u \cdot \nu = 3xyz exp(-t),                                x in \Gamma
%
% \Omega = unit sphere
%
% Exact solution:
%
% u(x,t) =   xyz exp(t)
% v(x,t) = 2 xyz exp(t)
%
%
% Note: if the solution is exponentially-in-time increasing instead of
% decreasing, it causes stability issues with the method. 
% The solution blows up. The stability condition probably becomes too
% restrictive.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [u, v, L2err_prod_rel] = solver_parabolic3d_bs_sphere(SurfaceElements, ElementsPlot, P, M, K, MS, KS, R, T, tau)

esol_u = @(P,t)   P(:,1).*P(:,2).*P(:,3)*exp(-t);
esol_v = @(P,t) 2*P(:,1).*P(:,2).*P(:,3)*exp(-t);

tic
NT = ceil(T/tau);
u =  esol_u(P,0);
v =  esol_v(R'*P,0);

MIT1 = M+tau*K;
MIT2 = MS+tau*KS;
perm1 = symamd(MIT1);
perm2 = symamd(MIT2);
[L1,U1] = lu(MIT1(perm1,perm1),'vector');
[L2,U2] = lu(MIT2(perm2,perm2),'vector');

for i=0:NT-1
    
   %F1 = M*(u - tau*esol_u(P, i*tau)) + 3*tau*R*MS*esol_u(R'*P,i*tau);
   %F1 = M*(u - tau*u) + 3*tau*R*MS*esol_u(R'*P,i*tau); % NO DIFFERENCE
   F1 = M*(u - tau*u) + tau*R*MS*(R'*u+v); % SMALL DIFFERENCE
   %F2 = MS*(v + tau*22*esol_u(R'*P,i*tau));
   F2 = MS*(v + tau*(7*v + 8*R'*u)); % BIG DIFFERENCE
   
   u(perm1) =  U1\(L1\F1(perm1));
   v(perm2) =  U2\(L2\F2(perm2));

end
toc

es_u = esol_u(P,T);
es_v = esol_v(R'*P,T);
err_u = u - es_u;
err_v = v - es_v;
L2err_prod = sqrt(err_u'*M*err_u + err_v'*MS*err_v);
normsol = sqrt(es_u'*M*es_u + es_v'*MS*es_v);
L2err_prod_rel = L2err_prod/normsol;
    
% Plotting Numerical Solution - Bulk Component u
figure
set(gcf, 'Color','white')
hold on
for i=1:length(ElementsPlot)
       plot(ElementsPlot(i), u(ElementsPlot(i).Pind)); 
end
view(3)
set(gca,'FontSize',18)
xlabel('x')
ylabel('y')
zlabel('z','rot',0)
axis equal
title('u')
colorbar
colormap parula

% Plotting Numerical Error - Bulk Component u
figure
set(gcf, 'Color','white')
hold on
for i=1:length(ElementsPlot)
       plot(ElementsPlot(i), esol_u(ElementsPlot(i).P,T) - u(ElementsPlot(i).Pind)); 
end
view(3)
set(gca,'FontSize',18)
xlabel('x')
ylabel('y')
zlabel('z','rot',0)
axis equal
title('u_{exact}-u')
colorbar
colormap parula

% Plotting Numerical Solution - Surface Component v
figure
set(gcf,'Color','white')
trisurf(SurfaceElements, P(:,1), P(:,2), P(:,3), R*v, 'FaceColor', 'interp')
view(3)
set(gca,'FontSize',18)
xlabel('x')
ylabel('y')
zlabel('z','rot',0)
axis equal
xlim([-0.5,1])
title('v')
colorbar
colormap parula
hold on
Ccirc = [-0.5, 0, 0];   % Center of circle 
Rcirc = sqrt(3)/2;      % Radius of circle 
teta = 0:0.01:2*pi;
xcirc = Ccirc(1) + zeros(size(teta)) ;
ycirc = Ccirc(2) + Rcirc*cos(teta);
zcirc = Ccirc(3) + Rcirc*sin(teta) ;
plot3(xcirc,ycirc,zcirc,'k','LineWidth',1.5)

% Plotting Numerical Error - Surface Component v
figure
set(gcf,'Color','white')
trisurf(SurfaceElements, P(:,1), P(:,2), P(:,3), esol_v(P,T) - R*v, 'FaceColor', 'interp')
view(3)
set(gca,'FontSize',18)
xlabel('x')
ylabel('y')
zlabel('z','rot',0)
axis equal
xlim([-0.5,1])
title('$v_{exact}-v$')
colorbar
colormap parula
hold on
Ccirc = [-0.5, 0, 0];   % Center of circle 
Rcirc = sqrt(3)/2;      % Radius of circle 
teta = 0:0.01:2*pi;
xcirc = Ccirc(1) + zeros(size(teta)) ;
ycirc = Ccirc(2) + Rcirc*cos(teta);
zcirc = Ccirc(3) + Rcirc*sin(teta) ;
plot3(xcirc,ycirc,zcirc,'k','LineWidth',1.5)


end