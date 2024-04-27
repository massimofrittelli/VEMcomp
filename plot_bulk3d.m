function [] = plot_bulk3d(ElementsPlot, u, titlestring)

for i=1:length(ElementsPlot)
   plot(ElementsPlot(i), u(ElementsPlot(i).Pind)); 
   hold on
end
axis equal tight
xlabel('x')
ylabel('y','rot',0)
zlabel('z','rot',0)
if nargin > 2
    title(titlestring,'Interpreter','latex')
end
set(gca,'FontSize',14)
colormap parula
colorbar

end