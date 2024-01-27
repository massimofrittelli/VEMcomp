function [v] = solver_elliptic_surf(D, alpha, g, P, MS, KS, R)

    LHS = D*KS + alpha*MS;
    RHS = MS*g(R'*P);
    v = LHS\RHS;

end