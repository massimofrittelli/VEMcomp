fun = @(P) P(:,1).^2 + P(:,2).^2 - 1;
[P, h, SquareElements, NonSquareElements] = generate_mesh_flat_domain(fun, 1, 20);
%fun = @(P) (-1/3+cos(P(:,1))).^2 + 0.15*P(:,1) + P(:,2).^2 -1;
%[P, h, SquareElements, NonSquareElements] = generate_mesh_flat_domain(fun, 3.5, 35);
%fun = @(P) max(abs(P(:,2)-1/2)-1/2, abs(P(:,1))-1-1/2*sin(2*pi*P(:,2)));
%[P, h, SquareElements, NonSquareElements] = generate_mesh_flat_domain(fun, 2, 41);
%fun = @(P) abs(P(:,1)).^(0.8) + abs(P(:,2)).^(0.8) - 1 + 0.1*sin(pi*P(:,1)) + 0.1*sin(2*pi*P(:,2));
%[P, h, SquareElements, NonSquareElements] = generate_mesh_flat_domain(fun, 1, 41);
Elements = [SquareElements; NonSquareElements];

figure
set(gcf,'color','white')
ii = 1:length(Elements);
for i=1:length(ii)
   plot(Elements(ii(i))); 
   hold on
end
colormap jet
caxis([-1.2,1.3]);
view(2)
axis equal tight
xlabel('x')
ylabel('y','rot',0)
set(gca,'FontSize',18)
%title('Jar-shaped mesh','FontSize',26)