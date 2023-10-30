% DESCRIPTION - GENERATES ILLUSTRATION OF STEP 1 OF MESH EXTRUSION ALGORITHM
% FOR THE BSVEM 3D ELLIPTIC PAPER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[P, h, boundarynode, Elements] = generate_mesh_sphere_bounding_box(6);

figure
set(gcf,'color','white')
hold on
for i=66:125
   plot(Elements(i)); 
end
view(3)
axis equal tight
xlabel('x')
ylabel('y')
zlabel('z','rot',0)
set(gca,'FontSize',18, 'Position', [0 0 1 1])

colormap jet
clim([-1.5,2]);

[x,y,z] = sphere(40);      %# Makes a 21-by-21 point sphere
r = 1;                 %# A radius value
surf(r.*z,r.*y,r.*x,0*x+1,'FaceColor','interp','EdgeColor','none','FaceAlpha',0.5);
lightangle(-40,20)
lighting gouraud
axis off