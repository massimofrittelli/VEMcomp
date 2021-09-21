% Script to test the stiffness matrix on the 3d sphere

load('sphere41.mat')
N = length(P);

% f1 = @(P) P(:,1).^2 + P(:,2).^2 + P(:,3).^2;
% integral1 = f1(P)'*K*f1(P);
% err1 = abs(16/5*pi-integral1);
% err_rel1 = err1*5/(16*pi);
% 
% f2 = @(P) P(:,1);
% integral2 = f2(P)'*K*f2(P);
% err2 = abs(4/3*pi-integral2);
% err_rel2 = err1*3/(4*pi);

rsquare = @(P) P(:,1).^2 + P(:,2).^2 + P(:,3).^2;
esol = @(P) (1- rsquare(P)).^2;
f = @(P) 4*(3-5*rsquare(P)) + esol(P);
tic
u = (K+M)\(M*f(P));
toc
es = esol(P);
err = u - es;
L2err = sqrt(err'*M*err);
normsol = sqrt(es'*M*es);
L2errel = L2err/normsol;

% perm = symamd(K+M);
% fperm = M(perm,perm)*f(P(perm,:));
% MKperm = K(perm,perm) + M(perm,perm);
% tic
% uperm = MKperm\fperm;
% u(perm,1) = uperm;
% toc
% es = esol(P);
% err = u - es;
% L2err = sqrt(err'*M*err);
% normsol = sqrt(es'*M*es);
% L2errel = L2err/normsol;

% Plotting Numerical solution
indsol = P(:,1) >= 0;
figure
set(gcf,'Color','white')
scatter3(P(indsol,1), P(indsol,2), P(indsol,3), 30, u(indsol,1),'filled')
view(3)
set(gca,'FontSize',18)
xlabel('x')
ylabel('y')
zlabel('z','rot',0)
axis equal
title('Numerical solution')
colorbar
