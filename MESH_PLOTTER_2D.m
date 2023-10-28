close all
clearvars
tol = 1e-10;
range = [-1.5,1.5; -1.5,1.5];
Nx = 20;
% fun = @(P) P(:,1).^2 + P(:,2).^2 - 1;
fun = @(P) P(:,1).^4 - 0.45*P(:,1).^2 + P(:,2).^2 - 1;

[P, h, BulkElements, SurfaceElements] = generate_mesh_2d(fun, range, Nx, tol);

figure
set(gcf,'color','white')
ii = 1:length(BulkElements);
for i=1:length(ii)
   plot(BulkElements(ii(i))); 
   hold on
end
for i=1:length(SurfaceElements)
    plot(P(SurfaceElements(i,:),1), P(SurfaceElements(i,:),2), 'b', 'LineWidth', 3)
end
colormap jet
clim([-1.5,2]);
view(2)
axis equal tight
xlabel('x')
ylabel('y','rot',0)
set(gca,'FontSize',18)