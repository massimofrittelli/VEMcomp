% Script to test the stiffness matrix on the 3d sphere

load('sphere21.mat')
N = length(P);

f1 = @(P) P(:,1).^2 + P(:,2).^2 + P(:,3).^2;
integral1 = f1(P)'*K*f1(P);
err1 = abs(16/5*pi-integral1);
err_rel1 = err1*5/(16*pi);

f2 = @(P) P(:,1);
integral2 = f2(P)'*K*f2(P);
err2 = abs(4/3*pi-integral2);
err_rel2 = err1*3/(4*pi);