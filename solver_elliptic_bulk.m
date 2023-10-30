function u = solver_elliptic_bulk(n,D,alpha,f,P,M,K,R,bcond)
    
    bulknodes = (1:length(M))';
    if strcmp(bcond, 'dir')
        boundarynodes = sum(R,2) == 1;
    else
        boundarynodes = [];
    end
    bulknodes(boundarynodes) = [];

    Mdir = M(bulknodes, bulknodes);
    Kdir = K(bulknodes, bulknodes);
    Nbulk = length(bulknodes);
    LHS = spalloc(n*Nbulk, n*Nbulk, n*nnz(Mdir+Kdir));
    for i=1:n
        LHS((i-1)*Nbulk+1:i*Nbulk, (i-1)*Nbulk+1:i*Nbulk) = D(i)*Kdir + alpha(i)*Mdir; %#ok
        rhs = M*f{i}(P);
        RHS((i-1)*Nbulk+1:i*Nbulk, :) = rhs(bulknodes);
    end
    uu = LHS\RHS;
    ubulk = zeros(Nbulk,n);
    for i=1:n
        ubulk(:,i) = uu((i-1)*Nbulk+1:i*Nbulk,1);
    end

    u = zeros(length(M),n);
    u(bulknodes,:) = ubulk;

end