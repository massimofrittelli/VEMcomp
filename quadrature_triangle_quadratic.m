function [XY,W] = quadrature_triangle_quadratic(v)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Symmetric Gaussian Quadrature of order N=2 for a triangular domain
% Reference: "HIGH DEGREE EFFICIENT SYMMETRICAL GAUSSIAN QUADRATURE RULES
% FOR THE TRIANGLE", D.A. Dunavant, 1985.


W = ones(3,1)'*abs(det([v, ones(3,1)]))/6;
XY = (3*eye(3) + ones(3))*v/6;
