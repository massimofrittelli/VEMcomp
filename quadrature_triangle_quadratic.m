function [XY,W] = quadrature_triangle_quadratic(v)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% quadrature_triangle.m - Gaussian Quadrature of order N=2 for a triangular domain
%
% This scripts computes the N^2=4 nodes and weights for a triangle with
% vertices given by the 3x2 vector v. The nodes are produced by collapsing
% the square to a triangle. 
%
% Sample Usage: 
%
% >>[X,Y,Wx,Wy]=triquad(8,[0 0; 0 2; 2 1])
% >>f=@(x,y) exp(x+y);
% >>Q=Wx'*feval(f,X,Y)*Wy;
%
% Reference:  J.N. Lyness, Ronald Cools, A Survey of Numerical Cubature
%             over Triangles (1994)
%             http://citeseer.ist.psu.edu/lyness94survey.html
%
% Written by: Greg von Winckel
% Contact: gregvw(at)math(dot)unm(dot)edu
% http://math.unm.edu/~gregvw
%
% Optimised for the case N=2 by Massimo Frittelli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%A = [1/3 1/15 1/35];
%B1 = 2/9;
%B = 6/25;
%ab = [1/3 2; 1/15 2/9; 1/35 6/25];
%y =  [-1;1]*sqrt(3);
%Lp = [-1;1]*3*sqrt(3)/2;
%s = sqrt(2/9);
%X = [1-sqrt(6); 1+sqrt(6)]/5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd = [ 1, 0, 0; -1, 0, 1; 0, 1,-1]*v; 

Y = [1-sqrt(3); 1+sqrt(3)]/2;
X = [6-sqrt(6); 6+sqrt(6)]/10;
[Y,X] = meshgrid(Y,X); 
Y = Y.*X;
X = X(:);
Y = Y(:);

Wy = [1;1]/2;
Wx = [9-sqrt(6); 9+sqrt(6)]/36;
Wx = abs(det(cd(2:3,:)))*Wx;
W = reshape(Wy*Wx',4,1);

XY = [cd(1,1)+cd(2,1)*X+cd(3,1)*Y, cd(1,2)+cd(2,2)*X+cd(3,2)*Y];
