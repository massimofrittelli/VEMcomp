% DESCRIPTION - solves parabolic surface problem on the 3D sphere, plots
% solution and computes error.

%
% v_t - LaplaceBeltrami(v) = 13u,   in \Gamma \times [0,T]
% v_0 = xyz,                        in \Gamma
%
% \Gamma = unit spherical surface
%
% Exact solution:
%
% v(x,t) = xyz exp(t)
%
%
% Note: if the solution is exponentially-in-time increasing instead of
% decreasing, it causes stability issues with the method. 
% The solution blows up. The stability condition probably becomes too
% restrictive.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [v] = solver_parabolic3d_surf_sphere(g, R, P, MS, KS, T, tau, v0)

tic
NT = ceil(T/tau);
v =  v0;

MIT2 = MS+tau*KS;
perm2 = symamd(MIT2);
[L2,U2] = lu(MIT2(perm2,perm2),'vector');

for i=0:NT-1
   
   F2 = MS*(v + tau*(0*v + g(R'*P, i*tau)));
   v(perm2) =  U2\(L2\F2(perm2));

end
toc

end