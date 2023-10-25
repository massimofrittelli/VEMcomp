function [XYZ,W] = quadrature_tetrahedron_quadratic(v)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Symmetric Gaussian Quadrature of order N=2 for a tetrahedral domain
% Reference: "SYMMETRIC GAUSSIAN QUADRATURE FORMULAE FOR TETRAHEDRONAL 
% REGIONS", Yu Jinyun, 1984.

W = ones(4,1)'*abs(det([v, ones(4,1)]))/24;
XYZ = ((0.5854101966249685-0.138196601125015)*eye(4) + 0.138196601125015*ones(4))*v;