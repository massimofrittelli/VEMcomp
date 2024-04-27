function [] = plot_surf2d(P,SurfaceElements,v,titlestring)

for i=1:length(SurfaceElements)
    patch(P(SurfaceElements(i,:),1), P(SurfaceElements(i,:),2), P(SurfaceElements(i,:),2), 'EdgeColor','interp','LineWidth',3,'FaceColor','none');
    hold on
end
axis equal tight
xlabel('x')
ylabel('y','rot',0)
if nargin > 2
    title(titlestring,'Interpreter','latex')
end
set(gca,'FontSize',14)
colormap parula
colorbar
caxis([min(v), max(v)])


end