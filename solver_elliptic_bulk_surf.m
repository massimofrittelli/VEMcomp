function [u,v] = solver_elliptic_bulk_surf(DOmega, DGamma, alpha, beta, f, g, P, M, MS, K, KS, R, gamma, delta)

N = length(M);
numsol = [DOmega*K + alpha*M - gamma*R*MS*R', -delta*R*MS;
          0*MS*R', DGamma*KS + beta*MS]\ ...
         [M*f(P); MS*g(R'*P)];
u = numsol(1:N,1);
v = numsol(N+1:end,1);

end