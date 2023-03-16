close all

Nx = 11;
xmin = -1;
xmax = 1;
ymin = -1;
ymax = 1;
zmin = 0;
zmax = 1;
range = [xmin, xmax; ymin, ymax; zmin, zmax];
tol = 1e-10;


R = 1;
r = 0.7;
L1 = 0.6;
L2 = 0.8;
b = (R^2-r^2)^2/(L1-L2);
a = -b*L2;
ff = @(x) (x <= L1).*(0*x-R^2) + (x>L1 & x<L2).*(-r^2-sqrt(a+b*x)) + (x>=L2).*(-r^2+0*x);
fun = @(P) ff(P(:,3)) + P(:,1).^2 + P(:,2).^2;

%fun = @(P) P(:,1).^2 + P(:,2).^2 + P(:,3).^2  -1;
%fun = @(P) (P(:,1).^2 + P(:,2).^2 + P(:,3).^2).^2 - (P(:,1).^2 + P(:,2).^2 + P(:,3).^2)+0.2;


[P, h, CubicElements, NonCubicElements, CubicElementsToPlot] = generate_mesh_3D_domain(fun, range, Nx,tol);
ElementsToPlot = [NonCubicElements(1:round(end/2)); CubicElementsToPlot];


figure
set(gcf,'color','white')
hold on
for i=1:length(ElementsToPlot)
   plot(ElementsToPlot(i)); 
end
view(3)
axis equal tight
xlabel('x')
ylabel('y')
zlabel('z','rot',0)
set(gca,'FontSize',18)

colormap jet
% colormap spring %(parabolic paper)
caxis([-1.5,2]);