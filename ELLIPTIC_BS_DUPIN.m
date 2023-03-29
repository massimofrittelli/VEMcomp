% DESCRIPTION - Solves an elliptic B-S toy model on the sphere, plots solution and
% computes errors.

close all

load('dupin_21.mat')
N = length(P); % Overall amount of nodes
NGamma = length(MS); % Amount of boundary nodes

esol_u = @(P) P(:,1).^2 + P(:,2).^2 + P(:,3).^2;
esol_v = @(P) 0*P(:,1) + 1;
f_u = @(P) P(:,1).^2 + P(:,2).^2 + P(:,3).^2 - 6;
f_v = @(P) P(:,1).^2 + P(:,2).^2 + P(:,3).^2 + 1;
nasty_term = @(P) 36*(9*(P(:,1).^2 + P(:,2).^2 + P(:,3).^2) + 261/100);
lev_fun_x = @(P) nasty_term(P).*P(:,1) -48*(6*P(:,1)-sqrt(39)/10);
lev_fun_y = @(P) nasty_term(P).*P(:,2) -2*3249/25*P(:,2);
lev_fun_z = @(P) nasty_term(P).*P(:,3);
flux_u = @(P) ((nasty_term(P).*P(:,1) -48*(6*P(:,1)-sqrt(39)/10))*2.*P(:,1) ...
    + (nasty_term(P).*P(:,2) -2*3249/25*P(:,2))*2.*P(:,2) ...
    + nasty_term(P)*2.*P(:,3).^2)./sqrt(lev_fun_x(P).^2 + lev_fun_y(P).^2 + lev_fun_z(P).^2) + 1;

RM = spalloc(N, NGamma, NGamma);
boundarynode = unique(EGamma(:));
RM(boundarynode, 1:NGamma) = speye(NGamma); % reduction matrix

tic
numsol = [K+M, RM*MS; MS*RM', KS+MS]\[M*f_u(P) + RM*MS*flux_u(P(boundarynode,:)); MS*f_v(P(boundarynode,:))];
u = numsol(1:N,1);
v = numsol(N+1:end,1);
toc
es_u = esol_u(P);
es_v = esol_v(P(boundarynode,:));
err_u = u - es_u;
err_v = v - es_v;
% L2err_u = sqrt(err_u'*M*err_u);
% L2err_v = sqrt(err_v'*MS*err_v);
% H1err_u = sqrt(err_u'*K*err_u);
% H1err_v = sqrt(err_v'*KS*err_v);

L2err_product_abs = sqrt(err_u'*M*err_u + err_v'*MS*err_v);
H1err_product_abs = sqrt(err_u'*(K+M)*err_u + err_v'*(KS+MS)*err_v);
L2norm_product = sqrt(es_u'*M*es_u + es_v'*MS*es_v);
H1norm_product = sqrt(es_u'*(K+M)*es_u + es_v'*(KS+MS)*es_v);
L2err_product_rel = L2err_product_abs/L2norm_product;
H1err_product_rel = H1err_product_abs/H1norm_product;


%xcut = range(1,1) + h/sqrt(3)*2;
Rcirc = sqrt(.1);      % Radius of circle 

% Plotting Numerical Solution - Bulk Component u
figure
set(gcf, 'Color','white')
hold on
for i=1:length(FacesToPlot)
       plot(FacesToPlot(i), u(FacesToPlot(i).Pind)); 
end
view(3)
set(gca,'FontSize',18)
xlabel('x')
ylabel('y')
zlabel('z','rot',0)
axis equal
title('$U$', 'interpreter', 'latex')
colorbar
colormap jet


% Plotting Numerical Solution - Surface Component v
% figure
% set(gcf,'Color','white')
% trisurf(EGamma, P(:,1), P(:,2), P(:,3), RM*v, 'EdgeColor', 'none', 'FaceColor', 'interp')
% view(3)
% set(gca,'FontSize',18)
% xlabel('x')
% ylabel('y')
% zlabel('z','rot',0)
% axis equal
% xlim([-0.41,1])
% title('$V$', 'interpreter', 'latex')
% colorbar
% colormap jet