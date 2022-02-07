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
% - Kbot,Mbot: stiffness and mass matrices on the bottom face
% - A: reduction matrix of the whole boundary
% - Abot: reduction matrix of the bottom face
% - Atop: reduction matrix of the entire cube without top face
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [P, h, K, M, Kbot, Mbot, KGamma, MGamma, A, Abot, Atop] = generate_mesh_cube(Nx)

% INITIALIZE VARIABLES
NGamma = 6*Nx^2 - 12*Nx + 8; % Number of nodes on the boundary
N = Nx^3; % Number of all nodes
hx = 1/(Nx-1); % Discretisation step along each dimension
h = hx*sqrt(3); % Meshsize
x = linspace(0,1,Nx);
P = zeros(N,3);
boundarynodecounter = 0;
xleftcounter = 0;
xrightcounter = 0;
yleftcounter = 0;
yrightcounter = 0;
zleftcounter = 0;
zrightcounter = 0;
boundarynode = zeros(NGamma,1); % Indexes of nodes that are on the boundary
xleftnode = zeros(Nx^2,1);
xrightnode = zeros(Nx^2,1);
yleftnode = zeros(Nx^2,1);
yrightnode = zeros(Nx^2,1);
zleftnode = zeros(Nx^2,1);
zrightnode = zeros(Nx^2,1);
A = spalloc(N, NGamma, NGamma); % reduction matrix of the whole boundary
Abot = spalloc(N,Nx^2,Nx^2); % Reduction matrix of the bottom face
Atop = spalloc(N,N-Nx^2,N-Nx^2); % Reduction matrix of cube without top face
Abot(1:Nx^2, 1:Nx^2) = speye(Nx^2);
Atop(1:N-Nx^2, 1:N-Nx^2) = speye(N-Nx^2);

% NODES OF REFERENCE CUBE
P1 = [0 0 0; 0 1 0; 1 1 0; 1 0 0]; % bottom face
P2 = [0 0 1; 0 1 1; 1 1 1; 1 0 1]; % top face
P3 = [0 0 0; 0 1 0; 0 1 1; 0 0 1]; % back face
P4 = [1 0 0; 1 1 0; 1 1 1; 1 0 1]; % front face
P5 = [0 0 0; 1 0 0; 1 0 1; 0 0 1]; % left face
P6 = [0 1 0; 1 1 0; 1 1 1; 0 1 1]; % right face

% LOCAL MATRICES OF REFERENCE CUBE
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

% LOCAL MATRICES OF REFERENCE SQUARE
KES = E1S.K;
MES = E1S.M;
KES([2 3], [2 3]) = KES([3 2], [3 2]); % reordering nodes
MES([2 3], [2 3]) = MES([3 2], [3 2]); % reordering nodes

% GLOBAL MATRICES IN THE BULK
K = spalloc(N,N,57*N);
M = spalloc(N,N,57*N);
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

% GLOBAL MATRICES ON THE BOTTOM FACE AND ON WHOLE SURFACE
Kbot = spalloc(Nx^2,Nx^2,9*N);
Mbot = spalloc(Nx^2,Nx^2,9*N);
KGamma = spalloc(N, N, 9*NGamma);
MGamma = spalloc(N, N, 9*NGamma);
% TODO: generate MGamma and KGamma

for j=0:Nx-2
    for k=0:Nx-2
            indexes = [Nx*j+k+1
                       Nx*j+k+2
                       Nx*(j+1)+k+1
                       Nx*(j+1)+k+2];
            Kbot(indexes, indexes) = Kbot(indexes, indexes) + KES; %#ok
            Mbot(indexes, indexes) = Mbot(indexes, indexes) + MES; %#ok
    end
end
Mbot = Mbot*hx^2;

% GENERATE NODES AND REDUCTION MATRIX OF WHOLE BOUNDARY
for i=0:Nx-1
    for j=0:Nx-1
        for k=0:Nx-1
            P(Nx^2*i+Nx*j+k+1,:) = [x(k+1) x(j+1) x(i+1)];
            if max(ismember([i j k], [0 Nx-1]))
                boundarynodecounter = boundarynodecounter+1;
                boundarynode(boundarynodecounter,1) = Nx^2*i+Nx*j+k+1;
                if i == 0
                    zleftcounter = zleftcounter + 1;
                    zleftnode(zleftcounter,1) = Nx^2*i+Nx*j+k+1;
                end
                if i == Nx-1
                    zrightcounter = zrightcounter + 1;
                    zrightnode(zrightcounter,1) = Nx^2*i+Nx*j+k+1;
                end
                if j == 0
                    yleftcounter = yleftcounter + 1;
                    yleftnode(yleftcounter,1) = Nx^2*i+Nx*j+k+1;
                end
                if j == Nx-1
                   yrightcounter = yrightcounter + 1; 
                   yrightnode(yrightcounter,1) = Nx^2*i+Nx*j+k+1;
                end
                if k == 0
                    xleftcounter = xleftcounter + 1;
                    xleftnode(xleftcounter,1) = Nx^2*i+Nx*j+k+1;
                end
                if k == Nx-1
                   xrightcounter = xrightcounter + 1; 
                   xrightnode(xrightcounter,1) = Nx^2*i+Nx*j+k+1;
                end
            end
        end
    end
end
A(boundarynode, 1:NGamma) = speye(NGamma);

MGamma(xleftnode, xleftnode) = MGamma(xleftnode, xleftnode) + Mbot;
MGamma(xrightnode, xrightnode) = MGamma(xrightnode, xrightnode) + Mbot;
MGamma(yleftnode, yleftnode) = MGamma(yleftnode, yleftnode) + Mbot;
MGamma(yrightnode, yrightnode) = MGamma(yrightnode, yrightnode) + Mbot;
MGamma(zleftnode, zleftnode) = MGamma(zleftnode, zleftnode) + Mbot;
MGamma(zrightnode, zrightnode) = MGamma(zrightnode, zrightnode) + Mbot;

KGamma(xleftnode, xleftnode) = KGamma(xleftnode, xleftnode) + Kbot;
KGamma(xrightnode, xrightnode) = KGamma(xrightnode, xrightnode) + Kbot;
KGamma(yleftnode, yleftnode) = KGamma(yleftnode, yleftnode) + Kbot;
KGamma(yrightnode, yrightnode) = KGamma(yrightnode, yrightnode) + Kbot;
KGamma(zleftnode, zleftnode) = KGamma(zleftnode, zleftnode) + Kbot;
KGamma(zrightnode, zrightnode) = KGamma(zrightnode, zrightnode) + Kbot;

MGamma = A'*MGamma*A;
KGamma = A'*KGamma*A;


end