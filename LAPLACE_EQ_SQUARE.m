% Script to test VEM 2d for the Laplace equation on the unit square

%   - Delta(u) + u = (2*pi^2+1)*cos(pi*x)*cos(pi*y)*cos(pi*z)

% Exact solution: u(x,y,z) = cos(pi*x)*cos(pi*y)*cos(pi*z)

clearvars

Nx = 41; % Number of gridpoints along each dimension
h = 1/(Nx-1); % Meshsize
x = linspace(0,1,Nx);
P = zeros(Nx^2,3);
for j=0:Nx-1
    for k=0:Nx-1
        P(Nx*j+k+1,:) = [x(j+1) x(k+1) 0];
    end
end
  
PE = [0 0 0; 0 1 0; 1 1 0; 1 0 0]; % bottom face

E = element2dsquare(PE);

KE = E.K;
ME = E.M;

K = spalloc(Nx^2,Nx^2,9*Nx^2);
M = spalloc(Nx^2,Nx^2,9*Nx^2);
    for j=0:Nx-2
        for k=0:Nx-2
            indexes = [Nx*j+k+1
                       Nx*j+k+2
                       Nx*(j+1)+k+2
                       Nx*(j+1)+k+1];
            K(indexes, indexes) = K(indexes, indexes) + KE; %#ok
            M(indexes, indexes) = M(indexes, indexes) + ME; %#ok
        end
    end
M = M*h^2;

f = @(x,y) (2*pi^2+1)*cos(pi*x).*cos(pi*y);
esol = @(x,y) cos(pi*x).*cos(pi*y);
es = esol(P(:,1),P(:,2));

% Solving linear system (K+M)u = Mf
u = (K+M)\(M*f(P(:,1),P(:,2)));

% Plotting exact solution
figure(1)
scatter(P(:,1), P(:,2), 30, esol(P(:,1), P(:,2)),'filled')
colorbar

% Plotting numerical solution
figure(2)
scatter(P(:,1), P(:,2), 30, u,'filled')
colorbar

% Computing L2 error
err = es - u;
l2err = sqrt(err'*M*err);
l2solnorm = sqrt(es'*M*es);
l2err_rel = l2err/l2solnorm;