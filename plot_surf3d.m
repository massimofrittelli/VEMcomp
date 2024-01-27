function [] = plot_surf3d(P,R,SurfaceElements,v,titlestring)

trisurf(SurfaceElements, P(:,1), P(:,2), P(:,3), R*v, 'FaceColor', 'interp')
axis equal tight
xlabel('x')
ylabel('y','rot',0)
zlabel('z','rot',0)
if nargin > 2
    title(titlestring,'Interpreter','latex')
end
set(gca,'FontSize',18)
colormap parula
colorbar
clim([min(v), max(v)])


end