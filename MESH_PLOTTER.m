% Plots the mesh generated as follows
% [P,h,K,M,radii,Elements] = generate_mesh_sphere(6);

figure
set(gcf,'color','white')
%ii = [73:207];
ii = 1:73;
hold on
for i=1:length(ii)
   plot(Elements(ii(i))); 
end
view(3)
axis equal
xlabel('x')
ylabel('y')
zlabel('z','rot',0)
title('Extruded cubic mesh on the sphere')
set(gca,'FontSize',18)