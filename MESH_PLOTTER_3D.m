close all

Nx = 21;
fun = @(P) 0.5*sin(2*pi*P(:,1)).^2 + 0.8*P(:,1).^2 + P(:,2).^2 + P(:,3).^2  -1;
%fun = @(P) (P(:,1).^2 + P(:,2).^2 + P(:,3).^2).^2 - (P(:,1).^2 + P(:,2).^2 + P(:,3).^2)+0.2;
xmax = 1.2;
tol = 1e-10;

[P, h, CubicElements, NonCubicElements] = generate_mesh_3D_domain(fun, xmax, Nx,tol);
ElementsToPlot = [CubicElements(1:ceil(end/2)); NonCubicElements(1:ceil(end/2))];


figure
set(gcf,'color','white')
hold on
for i=1:length(ElementsToPlot)
   plot(ElementsToPlot(i)); 
end
view(3)
%axis equal tight
xlabel('x')
ylabel('y')
zlabel('z','rot',0)
set(gca,'FontSize',18)

colormap jet
% colormap spring %(parabolic paper)
caxis([-1.5,2]);