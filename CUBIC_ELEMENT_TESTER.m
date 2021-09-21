clearvars

edge = 1; % edge length
shift = [0 0 0];

P1 = [0 0 0; 0 1 0; 1 1 0; 1 0 0]*edge + shift; % bottom face
P2 = [0 0 1; 0 1 1; 1 1 1; 1 0 1]*edge + shift; % top face
P3 = [0 0 0; 0 1 0; 0 1 1; 0 0 1]*edge + shift; % back face
P4 = [1 0 0; 1 1 0; 1 1 1; 1 0 1]*edge + shift; % front face
P5 = [0 0 0; 1 0 0; 1 0 1; 0 0 1]*edge + shift; % left face
P6 = [0 1 0; 1 1 0; 1 1 1; 0 1 1]*edge + shift; % right face

% NAIVE IMPLEMENTATION: ACTUALLY COMPUTE MATRICES
tic
E1 = element2d(P1, sum(P1,1)/4);
E2 = element2d(P2, sum(P2,1)/4);
E3 = element2d(P3, sum(P3,1)/4);
E4 = element2d(P4, sum(P4,1)/4);
E5 = element2d(P5, sum(P5,1)/4);
E6 = element2d(P6, sum(P6,1)/4);

P = unique([P1; P2; P3; P4; P5; P6],'rows');
E = element3d([E1;E2;E3;E4;E5;E6], P, sum(P,1)/8);
toc

% SMART IMPLEMENTATION: USE CLOSED FORM MATRICES FOR 2D
tic
E1S = element2dsquare(P1);
E2S = element2dsquare(P2);
E3S = element2dsquare(P3);
E4S = element2dsquare(P4);
E5S = element2dsquare(P5);
E6S = element2dsquare(P6);

PS = unique([P1; P2; P3; P4; P5; P6],'rows');
ES = element3dcube([E1S;E2S;E3S;E4S;E5S;E6S], PS);
toc

% CUT TEST
tic
E1D = element2ddummy(P1, true);
E2D = element2ddummy(P2, true);
E3D = element2ddummy(P3, true);
E4D = element2ddummy(P4, true);
E5D = element2ddummy(P5, true);
E6D = element2ddummy(P6, true);

PD = unique([P1; P2; P3; P4; P5; P6],'rows');
ED = element3ddummy(PD, [E1D;E2D;E3D;E4D;E5D;E6D], true);
EE = extrude(E2D);
toc