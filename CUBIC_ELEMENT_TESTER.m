clearvars

edge = 1; % edge length

P1 = [0 0 0; 0 1 0; 1 1 0; 1 0 0]*edge; % bottom face
P2 = [0 0 1; 0 1 1; 1 1 1; 1 0 1]*edge; % top face
P3 = [0 0 0; 0 1 0; 0 1 1; 0 0 1]*edge; % back face
P4 = [1 0 0; 1 1 0; 1 1 1; 1 0 1]*edge; % front face
P5 = [0 0 0; 1 0 0; 1 0 1; 0 0 1]*edge; % left face
P6 = [0 1 0; 1 1 0; 1 1 1; 0 1 1]*edge; % right face

% NAIVE IMPLEMENTATION: ACTUALLY COMPUTE MATRICES
tic
E1 = element2d(P1, sum(P1,1)/4);
E2 = element2d(P2, sum(P2,1)/4);
E3 = element2d(P3, sum(P3,1)/4);
E4 = element2d(P4, sum(P4,1)/4);
E5 = element2d(P5, sum(P5,1)/4);
E6 = element2d(P6, sum(P6,1)/4);
toc

P = unique([P1; P2; P3; P4; P5; P6],'rows');
E = element3d([E1;E2;E3;E4;E5;E6], P, [.5 .5 .5]);

% SMART IMPLEMENTATION: USE CLOSED FORM MATRICES FOR 2D
tic
E1S = element2dsquare(P1);
E2S = element2dsquare(P2);
E3S = element2dsquare(P3);
E4S = element2dsquare(P4);
E5S = element2dsquare(P5);
E6S = element2dsquare(P6);
toc

PS = unique([P1; P2; P3; P4; P5; P6],'rows');
ES = element3d([E1S;E2S;E3S;E4S;E5S;E6S], PS, [.5 .5 .5]);