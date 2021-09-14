% Script to test VEM 3d for the Laplace equation on the unit cube

%   - Delta(u) + u = (3*pi^2+1)*cos(pi*x)*cos(pi*y)*cos(pi*z)

% Exact solution: u(x,y,z) = cos(pi*x)*cos(pi*y)*cos(pi*z)

Nx = 41; % Number of gridpoints along each dimension
h = 1/(Nx-1); % Meshsize
x = linspace(0,1,Nx);
P = zeros(Nx^3,3);
for i=0:Nx-1
    for j=0:Nx-1
        for k=0:Nx-1
            P(Nx^2*i+Nx*j+k+1,:) = [x(i+1) x(j+1) x(k+1)];
        end
    end
end
  
P1 = [0 0 0; 0 1 0; 1 1 0; 1 0 0]; % bottom face
P2 = [0 0 1; 0 1 1; 1 1 1; 1 0 1]; % top face
P3 = [0 0 0; 0 1 0; 0 1 1; 0 0 1]; % back face
P4 = [1 0 0; 1 1 0; 1 1 1; 1 0 1]; % front face
P5 = [0 0 0; 1 0 0; 1 0 1; 0 0 1]; % left face
P6 = [0 1 0; 1 1 0; 1 1 1; 0 1 1]; % right face

E1 = element2d(P1, sum(P1,1)/4);
E2 = element2d(P2, sum(P2,1)/4);
E3 = element2d(P3, sum(P3,1)/4);
E4 = element2d(P4, sum(P4,1)/4);
E5 = element2d(P5, sum(P5,1)/4);
E6 = element2d(P6, sum(P6,1)/4);

PE = unique([P1; P2; P3; P4; P5; P6],'rows');
E = element3d([E1;E2;E3;E4;E5;E6], PE, sum(PE,1)/8);
KE = E.K;
ME = E.M;

K = spalloc(Nx^3,Nx^3,57*Nx^3);
M = spalloc(Nx^3,Nx^3,57*Nx^3);
for i=0:Nx-2
    for j=0:Nx-2
        for k=0:Nx-2
            indexes = [Nx^2*i+Nx*j+k+1
                       Nx^2*i+Nx*j+k+2
                       Nx^2*i+Nx*(j+1)+k+1
                       Nx^2*i+Nx*(j+1)+k+2
                       Nx^2*(i+1)+Nx*j+k+1
                       Nx^2*(i+1)+Nx*j+k+2
                       Nx^2*(i+1)+Nx*(j+1)+k+1
                       Nx^2*(i+1)+Nx*(j+1)+k+2];
            K(indexes, indexes) = K(indexes, indexes) + KE; %#ok
            M(indexes, indexes) = M(indexes, indexes) + ME; %#ok
        end
    end
end
K = K*h;
M = M*h^3;

f = @(x,y,z) (3*pi^2+1)*cos(pi*x).*cos(pi*y).*cos(pi*z);
esol = @(x,y,z) cos(pi*x).*cos(pi*y).*cos(pi*z);

% Solving linear system (K+M)u = Mf
tic
u = (K+M)\(M*f(P(:,1),P(:,2),P(:,3)));
es = esol(P(:,1),P(:,2),P(:,3));
toc

% Plotting Exact solution
figure(1)
scatter3(P(:,1), P(:,2), P(:,3), 30, esol(P(:,1), P(:,2), P(:,3)),'filled')
colorbar

% Plotting Numerical solution
figure(2)
scatter3(P(:,1), P(:,2), P(:,3), 30, u,'filled')
colorbar

% Computing L2 error
err = es - u;
l2err = sqrt(err'*M*err);