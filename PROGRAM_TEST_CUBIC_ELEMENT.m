% DESCRIPTION - Tests correctness of local matrices in 3d for the cube and
% extrusion process.

clearvars

edge = 1; % edge length
shift = [0 0 0]; % position of first vertex

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

Ksquare = [3/4           -1/4           -1/4           -1/4;     
      -1/4            3/4           -1/4           -1/4;     
      -1/4           -1/4            3/4           -1/4;     
      -1/4           -1/4           -1/4            3/4   ];
Msquare = [17/48          -3/16          13/48          -3/16;   
      -3/16          17/48          -3/16          13/48;    
      13/48          -3/16          17/48          -3/16;    
      -3/16          13/48          -3/16          17/48 ];
Csquare = [5/48           1/16           1/48           1/16;    
       1/16           5/48           1/16           1/48;    
       1/48           1/16           5/48           1/16;    
       1/16           1/48           1/16           5/48;    
];

max(max(abs(Ksquare - E1.K)))
max(max(abs(Csquare - E1.C)))
max(max(abs(Msquare - E1.M)))

P = unique([P1; P2; P3; P4; P5; P6],'rows');
E = element3d([E1;E2;E3;E4;E5;E6], P, sum(P,1)/8);

Kcube = [3 1 1 -1 1 -1 -1 -3;
1 3 -1 1 -1 1 -3 -1;
1 -1 3 1 -1 -3 1 -1;
-1 1 1 3 -3 -1 -1 1;
1 -1 -1 -3 3 1 1 -1;
-1 1 -3 -1 1 3 -1 1
-1 -3 1 -1 1 -1 3 1
-3 -1 -1 1 -1 1 1 3]/16 ...
+ [2 -1 -1 0 -1 0 0 1
-1 2 0 -1 0 -1 1 0
-1 0 2 -1 0 1 -1 0
0 -1 -1 2 1 0 0 -1
-1 0 0 1 2 -1 -1 0
0 -1 1 0 -1 2 0 -1
0 1 -1 0 -1 0 2 -1
1 0 0 -1 0 -1 -1 2]*sqrt(3)/4;

Mcube = [17/32         -11/48         -11/48           1/96         -11/48           1/96           1/96           1/4;     
     -11/48          17/32           1/96         -11/48           1/96         -11/48           1/4            1/96;    
     -11/48           1/96          17/32         -11/48           1/96           1/4          -11/48           1/96;    
       1/96         -11/48         -11/48          17/32           1/4            1/96           1/96         -11/48;    
     -11/48           1/96           1/96           1/4           17/32         -11/48         -11/48           1/96;    
       1/96         -11/48           1/4            1/96         -11/48          17/32           1/96         -11/48;    
       1/96           1/4          -11/48           1/96         -11/48           1/96          17/32         -11/48;    
       1/4            1/96           1/96         -11/48           1/96         -11/48         -11/48          17/32   
];

Ccube = [1/32           1/48           1/48           1/96           1/48           1/96           1/96           0;     
       1/48           1/32           1/96           1/48           1/96           1/48           0              1/96;   
       1/48           1/96           1/32           1/48           1/96           0              1/48           1/96;   
       1/96           1/48           1/48           1/32           0              1/96           1/96           1/48;   
       1/48           1/96           1/96           0              1/32           1/48           1/48           1/96;   
       1/96           1/48           0              1/96           1/48           1/32           1/96           1/48;   
       1/96           0              1/48           1/96           1/48           1/96           1/32           1/48;   
       0              1/96           1/96           1/48           1/96           1/48           1/48           1/32    
];

max(max(abs(Kcube - E.K)))
max(max(abs(Ccube - E.C)))
max(max(abs(Mcube - E.M)))

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
PD = unique([P1; P2; P3; P4; P5; P6],'rows');

E1D = element2d_dummy(P1, true, false, find(ismember(PD, P1,'rows')));
E2D = element2d_dummy(P2, true, false, find(ismember(PD, P2,'rows'))); 
E3D = element2d_dummy(P3, true, false, find(ismember(PD, P3,'rows'))); 
E4D = element2d_dummy(P4, true, false, find(ismember(PD, P4,'rows')));
E5D = element2d_dummy(P5, true, false, find(ismember(PD, P5,'rows')));
E6D = element2d_dummy(P6, true, false, find(ismember(PD, P6,'rows')));

ED = element3d_dummy(PD, [E1D;E2D;E3D;E4D;E5D;E6D], true, 1:8);
EE = extrude(E2D, 8);
toc