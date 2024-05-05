% This scripts tests the correctness of the numerically computed 
% VEM local matrices for the unit cube

% Create the unit cube E and compute its local VEM matrices
P1 = [0 0 0; 0 1 0; 1 1 0; 1 0 0]; % bottom face
E1 = element2d(P1,  false, false, [], sum(P1,1)/4);
P2 = [0 0 1; 0 1 1; 1 1 1; 1 0 1]; % top face
E2 = element2d(P2,  false, false, [],  sum(P2,1)/4);
P3 = [0 0 0; 0 1 0; 0 1 1; 0 0 1]; % back face
E3 = element2d(P3,  false, false, [],  sum(P3,1)/4);
P4 = [1 0 0; 1 1 0; 1 1 1; 1 0 1]; % front face
E4 = element2d(P4,  false, false, [],  sum(P4,1)/4);
P5 = [0 0 0; 1 0 0; 1 0 1; 0 0 1]; % left face
E5 = element2d(P5,  false, false, [],  sum(P5,1)/4);
P6 = [0 1 0; 1 1 0; 1 1 1; 0 1 1]; % right face
E6 = element2d(P6,  false, false, [],  sum(P6,1)/4);
P = unique([P1; P2; P3; P4; P5; P6],'rows');
E = getLocalMatrices(element3d(P, [E1;E2;E3;E4;E5;E6], false, [], sum(P,1)/8));
% Closed-form local VEM matrices
K = [3 1 1 -1 1 -1 -1 -3;
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
M = [17/32         -11/48         -11/48           1/96         -11/48           1/96           1/96           1/4;     
     -11/48          17/32           1/96         -11/48           1/96         -11/48           1/4            1/96;    
     -11/48           1/96          17/32         -11/48           1/96           1/4          -11/48           1/96;    
       1/96         -11/48         -11/48          17/32           1/4            1/96           1/96         -11/48;    
     -11/48           1/96           1/96           1/4           17/32         -11/48         -11/48           1/96;    
       1/96         -11/48           1/4            1/96         -11/48          17/32           1/96         -11/48;    
       1/96           1/4          -11/48           1/96         -11/48           1/96          17/32         -11/48;    
       1/4            1/96           1/96         -11/48           1/96         -11/48         -11/48          17/32   
];
C = [1/32           1/48           1/48           1/96           1/48           1/96           1/96           0;     
       1/48           1/32           1/96           1/48           1/96           1/48           0              1/96;   
       1/48           1/96           1/32           1/48           1/96           0              1/48           1/96;   
       1/96           1/48           1/48           1/32           0              1/96           1/96           1/48;   
       1/48           1/96           1/96           0              1/32           1/48           1/48           1/96;   
       1/96           1/48           0              1/96           1/48           1/32           1/96           1/48;   
       1/96           0              1/48           1/96           1/48           1/96           1/32           1/48;   
       0              1/96           1/96           1/48           1/96           1/48           1/48           1/32    
];
% Evaluate maximum error of numerically computed local matrices
norm(K(:) - E.K(:),'inf')
norm(C(:) - E.C(:),'inf')
norm(M(:) - E.M(:),'inf')