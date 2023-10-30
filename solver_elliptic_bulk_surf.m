function [u] = solver_elliptic_bulk_surf(n, D, alpha, f, P, M, K)

    Nbulk = length(M);
    LHS = spalloc(n*Nbulk, n*Nbulk, n*nnz(M+K));
    for i=1:n
        LHS((i-1)*Nbulk+1:i*Nbulk, (i-1)*Nbulk+1:i*Nbulk) = D(i)*K + alpha(i)*M; %#ok
        RHS((i-1)*Nbulk+1:i*Nbulk, :) = M*f{i}(P);
    end
    uu = LHS\RHS;
    u = zeros(Nbulk,n);
    for i=1:n
        u(:,i) = uu((i-1)*Nbulk+1:i*Nbulk,1);
    end

end