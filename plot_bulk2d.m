function [] = plot_bulk2d(BulkElements, u, titlestring)

for i=1:length(BulkElements)
   plot(BulkElements(i), u(BulkElements(i).Pind),'none'); 
   hold on
end
view(2)
axis equal tight
xlabel('x')
ylabel('y','rot',0)
if nargin > 2
    title(titlestring, 'Interpreter', 'latex')
end
set(gca,'FontSize',14)
colormap parula
colorbar


end