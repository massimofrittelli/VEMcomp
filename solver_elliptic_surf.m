function [u] = solver_elliptic_surf(n, D, alpha, f, P, MS, KS, R)

    Nsurf = length(MS);
    LHS = spalloc(n*Nsurf, n*Nsurf, n*nnz(MS+KS));
    for i=1:n
        LHS((i-1)*Nsurf+1:i*Nsurf, (i-1)*Nsurf+1:i*Nsurf) = D(i)*KS + alpha(i)*MS; %#ok
        RHS((i-1)*Nsurf+1:i*Nsurf, :) = MS*R'*f{i}(P);
    end
    uu = LHS\RHS;
    u = zeros(Nsurf,n);
    for i=1:n
        u(:,i) = uu((i-1)*Nsurf+1:i*Nsurf,1);
    end

end