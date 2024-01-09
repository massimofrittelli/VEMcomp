function [] = plot_bulk_2d(BulkElements, u, titlestring)

for i=1:length(BulkElements)
   plot(BulkElements(i), u(BulkElements(i).Pind)); 
   hold on
end
view(2)
axis equal tight
xlabel('x')
ylabel('y','rot',0)
if nargin > 2
    title(titlestring, 'Interpreter', 'latex')
end
set(gca,'FontSize',18)
colormap jet
colorbar


end