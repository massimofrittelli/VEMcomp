function [X,Y,Z,W] = quadrature_tetrahedron_quadratic(vert)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Gaussian Quadrature for a tetrahedron of order N = 2
%
% Construct Gauss points and weights for a tetrahedral domain with vertices
% specified by the 4x3 matrix vert. Where each row contains the (x,y,z) for 
% a vertex.
%
% Sample usage: 
%
% Suppose f(x,y,z)=x^2+y^2+z^2. Then let's integrate this over a regular
% tetrahedron.
%
% >>vert=[1/sqrt(3) 0 0; -sqrt(3)/6,1/2,0;-sqrt(3)/6,-1/2,0;0 0 sqrt(6)/3];
% >>[X,Y,Z,W]=tetraquad(4,vert);
% >>F=X.^2+Y.^2+Z.^2;
% >>Q=W'*F;
%
% Written by: Greg von Winckel  
% Contact: gregvw(at)math(dot)unm(dot)edu
% http://math.unm.edu/~gregvw
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q1 = [10-sqrt(10), 10+sqrt(10)]/15;
w1 = [8-sqrt(10); 8+sqrt(10)]/48;
q2 = [6-sqrt(6); 6+sqrt(6)]/10; 
w2 = [9-sqrt(6); 9+sqrt(6)]/36;
q3 = [-sqrt(1/3) + 1; sqrt(1/3) + 1]/2;
w3 = [1;1]/2;
[q1,q2,q3]=meshgrid(q1,q2,q3);
q1=q1(:);   q2=q2(:);       q3=q3(:);
x=1-q1;     y=(1-q2).*q1;   z=q1.*q2.*q3;
w=reshape(reshape(w2*w1',4,1)*w3',8,1);
c=[1 0 0 0;-1 1 0 0;-1 0 1 0;-1 0 0 1]*vert;
W=abs(det(c(2:4,:)))*w;
% Change of coordinates 
XYZ=[ones(8,1) x y z]*c; X=XYZ(:,1); Y=XYZ(:,2); Z=XYZ(:,3);
 
% function [x,w]=rquad2
% k2 = 4;
% k1 = 3;
% A=[1/2 1/6 1/12];
% B1=3/20; 
% B=64/315;
% ab = [1/2 8/3; 1/6 3/20; 1/12 64/315];
% s = sqrt(3/20);
% [V,~]=eig([1/2 sqrt(3/20); sqrt(3/20) 1/6]);
% X = [5-sqrt(40), 5+sqrt(40)]/15;
% Grid points
% x = [10-sqrt(10), 10+sqrt(10)]/15;
% Quadrature weights
% w = [8-sqrt(10); 8+sqrt(10)]/48;

% function [x,w]=rquad1
% k1=2; 
% k2=3; 
% A = [1/3 1/15 1/35];
% B1 = 2/9; 
% B = 6/25;
% ab = [1/3 2; 1/15 2/9; 1/35 6/25]; 
% s = sqrt(2/9);
% [V,~] = eig([1/3 sqrt(2/9); sqrt(2/9) 1/15]);
% Grid points
% x = [6-sqrt(6); 6+sqrt(6)]/10; 
% Quadrature weights
% w = [9-sqrt(6); 9+sqrt(6)]/36;

% function [x,w]=rquad0
% k1=k+1; 
% k2=k+2; 
% A=[0 0 0];
% B1 = 1/3;
% B = 4/15;
% ab=[0 2; 0 1/3; 0 4/15]; 
% s=sqrt(1/3);
% [V,~] = eig([0 sqrt(1/3); sqrt(1/3) 0]);
% X = [-sqrt(1/3); sqrt(1/3)];
% Grid points
% x = [-sqrt(1/3) + 1; sqrt(1/3) + 1]/2; 
% Quadrature weights
% w = [1;1]/2;