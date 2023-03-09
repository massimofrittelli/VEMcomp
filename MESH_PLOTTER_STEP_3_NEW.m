% DESCRIPTION - GENERATES ILLUSTRATION OF STEP 3 OF MESH EXTRUSION ALGORITHM
% FOR THE BSVEM 3D ELLIPTIC PAPER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%[P, h, K, M, KS, MS, boundarynode, EGamma, Elements] = plot_mesh_step_3(8);
[P, h, K, M, KS, MS, boundarynode, EGamma, Elements] = plot_mesh_step_3_new(13);
%[P, h, K, M, KS, MS, CMS, boundarynode, EGamma, Elements, EGammaCut] = generate_mesh_sphere_new(13);
figure
set(gcf,'color','white')
%ii = 101:207;
%ii = 40:73;
ii = 1:66;
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