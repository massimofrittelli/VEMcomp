% CREATES FIGURE THAT ILLUSTRATES THE DIB SQUARE DOMAIN WITH ITS BCs

close all

figure('Renderer', 'zbuffer', 'color', 'white')
hold on

h1 = fill3([0 0 1 1],[0 1 1 0],[1 1 1 1],'b');
h2 = fill3([0 0 1 1],[0 0 0 0],[0 1 1 0],'g');
h3 = fill3([0 0 1 1],[1 1 1 1],[0 1 1 0],'g');
h4 = fill3([0 0 0 0],[0 0 1 1],[0 1 1 0],'g');
h5 = fill3([1 1 1 1],[0 0 1 1],[0 1 1 0],'g');
h6 = fill3([0 0 1 1],[0 1 1 0],[0 0 0 0],'r');

h1.FaceAlpha = 0.5;
h2.FaceAlpha = 0.3;
h3.FaceAlpha = 0.3;
h4.FaceAlpha = 0.3;
h5.FaceAlpha = 0.3;
h6.FaceAlpha = 1;

plot3([0 0 0 0 0],[0 0 1 1 0],[0 1 1 0 0],'k', 'LineWidth',3)
plot3([0 0 1 1 0],[0 1 1 0 0],[1 1 1 1 1],'k', 'LineWidth',3)
plot3([0 0 1 1 0],[0 0 0 0 0],[0 1 1 0 0],'k', 'LineWidth',3)

legend([h1, h2, h6], {'$\Gamma_D$', '$\Gamma_L$','$\Gamma$'}, 'Interpreter', 'latex')

view(3)

set(gca, 'FontSize', 22)
axis equal
axis off