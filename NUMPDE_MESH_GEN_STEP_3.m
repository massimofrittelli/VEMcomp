% DESCRIPTION - GENERATES ILLUSTRATION OF STEP 3 OF MESH EXTRUSION ALGORITHM
% FOR THE BSVEM 3D ELLIPTIC PAPER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NEW ALGORITHM: EXTRUSION WITH QUADRILATERAL LATERAL FACES (DOES NOT LEAVE
% HOLES, DOES NOT HAVE DUPLICATE NODES, AND HAS TWO NODE EXTRUSION OPTIONS)
% - option 1: radial projection. This mahes mesh more uniform, but extruded
% lateral faces are not guaranteed to be flat
% - option 2: orthant-wise constant extrusion direction. This makes mesh
% less uniform, but all internal faces are guaranteed to be flat.
Nx = 6;
xmin = -1;
xmax = 1;
ymin = -1;
ymax = 1;
zmin = -1;
zmax = 1;
range = [xmin, xmax; ymin, ymax; zmin, zmax];
tol = 1e-10;
fun = @(P) P(:,1).^2 + P(:,2).^2 + P(:,3).^2  -1;
[P, h, BulkElements, SurfaceElements, ElementsPlot] = generate_mesh3d(fun, range, Nx, tol, -0.2);

figure
set(gcf,'color','white')
hold on
for i=1:length(ElementsPlot)
    if i == 13
        plot(ElementsPlot(i), 1.5); 
    else
        plot(ElementsPlot(i)); 
    end
end
view(3)
axis equal tight
xlabel('x')
ylabel('y')
zlabel('z','rot',0)
set(gca,'FontSize',18, 'Position', [0 0 1 1])

colormap jet
% colormap spring %(parabolic paper)
clim([-1.5,2]);
lightangle(-40,20)
lighting gouraud
axis off