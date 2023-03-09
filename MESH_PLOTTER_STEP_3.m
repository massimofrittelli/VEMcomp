% DESCRIPTION - GENERATES ILLUSTRATION OF STEP 3 OF MESH EXTRUSION ALGORITHM
% FOR THE BSVEM 3D ELLIPTIC PAPER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% OLD ALGORITHM: EXTRUSION WITH TRIANGULAR LATERAL FACES (LEAVES HOLES AND
% SUFFERS FROM DUPLICATED NODES)
%[P, h, K, M, KS, MS, boundarynode, EGamma, Elements] = plot_mesh_step_3(8);

% NEW ALGORITHM: EXTRUSION WITH QUADRILATERAL LATERAL FACES (DOES NOT LEAVE
% HOLES, DOES NOT HAVE DUPLICATE NODES, AND HAS TWO NODE EXTRUSION OPTIONS)
% - option 1: radial projection. This mahes mesh more uniform, but extruded
% lateral faces are not guaranteed to be flat
% - option 2: orthant-wise constant extrusion direction. This makes mesh
% less uniform, but all internal faces are guaranteed to be flat.
[P, h, boundarynode, EGamma, Elements] = plot_mesh_step_3_new(5);

figure
set(gcf,'color','white')
ii = 1:16;
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
% colormap spring %(parabolic paper)
caxis([-1.2,1.3]);

%lightangle(-40,20)
%lighting gouraud
axis off