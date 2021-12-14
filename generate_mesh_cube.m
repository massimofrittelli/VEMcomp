%
% Generates polyhedral mesh and matrices on the unit "DIB cube"
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INPUTS:
%
% - Nx: number of nodes along each direction of the cube
%
% OUTPUTS:
%
% - P: array of nodes
% - h: meshsize
% - K,M: stiffness and mass matrices in the bulk
% - KS,MS: stiffness and mass matrices on the bottom face
% - Abot: reduction matrix of the bottom bace
% - Atop: reduction matrix of the entire cube without top face
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [P, h, K, M, KS, MS, Abot, Atop] = generate_mesh_cube(Nx)

hx = 1/(Nx-1); % Discretisation step along each dimension
h = hx*sqrt(3); % Meshsize
x = linspace(0,1,Nx);
P = zeros(Nx^3,3);
Abot = spalloc(Nx^3,Nx^2,Nx^2);
Atop = spalloc(Nx^3,Nx^2-Nx^2,Nx^3-Nx^2);
Abot(1:Nx^2, 1:Nx^2) = speye(Nx^2);
Atop(1:Nx^3-Nx^2, 1:Nx^3-Nx^2) = speye(Nx^3-Nx^2);

for i=0:Nx-1
    for j=0:Nx-1
        for k=0:Nx-1
            P(Nx^2*i+Nx*j+k+1,:) = [x(k+1) x(j+1) x(i+1)];
        end
    end
end
  
P1 = [0 0 0; 0 1 0; 1 1 0; 1 0 0]; % bottom face
P2 = [0 0 1; 0 1 1; 1 1 1; 1 0 1]; % top face
P3 = [0 0 0; 0 1 0; 0 1 1; 0 0 1]; % back face
P4 = [1 0 0; 1 1 0; 1 1 1; 1 0 1]; % front face
P5 = [0 0 0; 1 0 0; 1 0 1; 0 0 1]; % left face
P6 = [0 1 0; 1 1 0; 1 1 1; 0 1 1]; % right face

E1S = element2dsquare(P1);
E2S = element2dsquare(P2);
E3S = element2dsquare(P3);
E4S = element2dsquare(P4);
E5S = element2dsquare(P5);
E6S = element2dsquare(P6);
PS = [0     0     0;
      1     0     0;
      0     1     0;
      1     1     0;
      0     0     1;
      1     0     1;
      0     1     1;
      1     1     1];
ES = element3dcube([E1S;E2S;E3S;E4S;E5S;E6S], PS);

KE = ES.K;
ME = ES.M;

KES = E1S.K;
MES = E1S.M;
KES([2 3], [2 3]) = KES([3 2], [3 2]);
MES([2 3], [2 3]) = MES([3 2], [3 2]);

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
K = K*hx;
M = M*hx^3;

KS = spalloc(Nx^2,Nx^2,9*Nx^3);
MS = spalloc(Nx^2,Nx^2,9*Nx^3);
for j=0:Nx-2
    for k=0:Nx-2
            indexes = [Nx*j+k+1
                       Nx*j+k+2
                       Nx*(j+1)+k+1
                       Nx*(j+1)+k+2];
            KS(indexes, indexes) = KS(indexes, indexes) + KES; %#ok
            MS(indexes, indexes) = MS(indexes, indexes) + MES; %#ok
    end
end
MS = MS*hx^2;

end