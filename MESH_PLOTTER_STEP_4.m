% DESCRIPTION - GENERATES ILLUSTRATION OF STEP 4 OF MESH EXTRUSION ALGORITHM
% FOR THE BSVEM 3D ELLIPTIC PAPER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NEW ALGORITHM: EXTRUSION WITH QUADRILATERAL LATERAL FACES (DOES NOT LEAVE
% HOLES, DOES NOT HAVE DUPLICATE NODES, AND HAS TWO NODE EXTRUSION OPTIONS)
% - option 1: radial projection. This mahes mesh more uniform, but extruded
% lateral faces are not guaranteed to be flat
% - option 2: orthant-wise constant extrusion direction. This makes mesh
% less uniform, but all internal faces are guaranteed to be flat.
[P, h, boundarynode, EGamma, Elements] = generate_mesh_sphere(9);

figure
set(gcf,'color','white')
hold on
i = 2;
plot(Elements(i)); 
view(3)
axis equal tight
xlabel('x')
ylabel('y')
zlabel('z','rot',0)
set(gca,'FontSize',18, 'Position', [0 0 1 1])

colormap jet
caxis([-1.2,1.3]);

axis off