close all

Nx = 11;
tol = 1e-10;

% BATTERY-SHAPED DOMAIN
% R = 1;
% r = 0.7;
% L1 = 0.6;
% L2 = 0.8;
% b = (R^2-r^2)^2/(L1-L2);
% a = -b*L2;
% ff = @(x) (x <= L1).*(0*x-R^2) + (x>L1 & x<L2).*(-r^2-sqrt(a+b*x)) + (x>=L2).*(-r^2+0*x);
% fun = @(P) ff(P(:,3)) + P(:,1).^2 + P(:,2).^2;

% UNIT BALL
% xmin = -1.04;
% xmax = 1.04;
% ymin = -1.04;
% ymax = 1.04;
% zmin = -1.04;
% zmax = 1.04;
% fun = @(P) P(:,1).^2 + P(:,2).^2 + P(:,3).^2  -1;

% HOLLOW BALL
% fun = @(P) (P(:,1).^2 + P(:,2).^2 + P(:,3).^2).^2 - (P(:,1).^2 + P(:,2).^2 + P(:,3).^2)+0.2;

% ALIEN SHAPED DOMAIN
% xmin = -1.6;
% xmax = 1.6;
% ymin = -1.2;
% ymax = 1.5;
% zmin = -1.2;
% zmax = 1.2;
% fun = @(P) P(:,1).^2 + .4*cos(2*pi*P(:,1)) + P(:,2).^2 + .4*sin(pi*P(:,2)) + P(:,3).^2 + .6*cos(3*pi*P(:,3)) - 1;

% DUMBBELL SHAPED DOMAIN
% xmin = -1.1;
% xmax = 1.1;
% ymin = -0.7;
% ymax = 0.7;
% zmin = -0.7;
% zmax = 0.7;
% fun = @(P) P(:,1).^4 - P(:,1).^2 + P(:,2).^2 + P(:,3).^2 - 0.1;

% DUPIN RING CYCLIDE
% xmin = -1.25;
% xmax = 1;
% ymin = -1.05;
% ymax = 1.05;
% zmin = -0.6;
% zmax = 0.6;
% fun = @(P) (9*(P(:,1).^2 + P(:,2).^2 + P(:,3).^2) + 261/100).^2 -4*(6*P(:,1)-sqrt(39)/10).^2 -3249/25*P(:,2).^2;

% TORUS
xmin = -1;
xmax = 1;
ymin = -1;
ymax = 1;
zmin = -0.31;
zmax = 0.31;
fun = @(P) (sqrt(P(:,1).^2 + P(:,2).^2) - 7/10).^2 + P(:,3).^2 - 9/100;

range = [xmin, xmax; ymin, ymax; zmin, zmax];
xcut = -.5;
[P, h, BulkElements, SurfaceElements, ElementsPlot] = generate_mesh_3d(fun, range, Nx, tol, xcut);

figure
set(gcf,'color','white')
hold on
for i=1:length(ElementsPlot)
    if ElementsPlot(i).is_boundary
        plot(ElementsPlot(i), -0.5); 
    else
        plot(ElementsPlot(i));
    end
end
view(3)
axis equal tight
xlabel('x')
ylabel('y','rot',0)
zlabel('z','rot',0)
set(gca,'FontSize',18)
colormap jet
% colormap spring %(parabolic paper)
clim([-1.5,2]);

% x = linspace(xmin, xmax, Nx);
% y = linspace(ymin, ymax, Nx);
% z = linspace(zmin, zmax, Nx);
% [X,Y,Z] = meshgrid(x,y,z);
% fun = @(X,Y,Z) (9*(X.^2 + Y.^2 + Z.^2) + 261/100).^2 -4*(6*X-sqrt(39)/10).^2 -3249/25*Y.^2;
% figure
% set(gcf,'color','white')
% isosurface(X, Y, Z, fun(X,Y,Z), 1);
% axis equal tight
% view(3)
% xlabel('x')
% ylabel('y')
% zlabel('z','rot',0)
% set(gca,'FontSize',18)
% colormap jet
% % colormap spring %(parabolic paper)
% caxis([-1.5,2]);
% %lightangle(-40,20)
% % %lighting gouraud
% 
