% DESCRIPTION - solves parabolic surface problem on the 3D sphere, plots
% solution and computes error.

%
% v_t - LaplaceBeltrami(v) = 13u,   in \Gamma \times [0,T]
% v_0 = xyz,                        in \Gamma
%
% \Gamma = unit spherical surface
%
% Exact solution:
%
% v(x,t) = xyz exp(t)
%
%
% Note: if the solution is exponentially-in-time increasing instead of
% decreasing, it causes stability issues with the method. 
% The solution blows up. The stability condition probably becomes too
% restrictive.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [v, L2err_rel] = solver_parabolic3d_surf_sphere(SurfaceElements, R, P, MS, KS, T, tau)

esol_v = @(P,t) P(:,1).*P(:,2).*P(:,3)*exp(t);

tic
NT = ceil(T/tau);
v =  esol_v(R'*P,0);

MIT2 = MS+tau*KS;
perm2 = symamd(MIT2);
[L2,U2] = lu(MIT2(perm2,perm2),'vector');

for i=0:NT-1
   
   F2 = MS*(v + tau*13*v);
   v(perm2) =  U2\(L2\F2(perm2));

end
toc

es_v = esol_v(R'*P,T);
err_v = v - es_v;
L2err_abs = sqrt(err_v'*MS*err_v);
normsol = sqrt(es_v'*MS*es_v);
L2err_rel = L2err_abs/normsol;
    

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
title('$v_{\text{exact}}-v$')
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