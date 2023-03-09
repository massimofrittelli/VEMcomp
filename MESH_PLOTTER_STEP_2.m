% DESCRIPTION - GENERATES ILLUSTRATION OF STEP 2 OF MESH EXTRUSION ALGORITHM
% FOR THE BSVEM 3D ELLIPTIC PAPER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%[P, h, K, M, KS, MS, boundarynode, EGamma, Elements] = plot_mesh_step_2(8);
[P, h, boundarynode, EGamma, Elements] = plot_mesh_step_2(6);

figure
set(gcf,'color','white')
%ii = 18:81;
ii = 1:19;
hold on
for i=1:length(ii)
   plot(Elements(ii(i))); 
end
view(3)
axis equal tight
xlabel('x')
ylabel('y')
zlabel('z','rot',0)
set(gca,'FontSize',18, 'Position', [0 0 1 1])


colormap jet
caxis([-1.2,1.3]);

[x,y,z] = sphere(40);      %# Makes a 21-by-21 point sphere
%x = x(21:end,:);       %# Keep top 11 x points
%y = y(21:end,:);       %# Keep top 11 y points
%z = z(21:end,:);       %# Keep top 11 z points
r = 1;                 %# A radius value
surf(r.*z,r.*y,r.*x,0*x+1,'FaceColor','interp','EdgeColor','none','FaceAlpha',0.5);

lightangle(-40,20)
lighting gouraud
axis off