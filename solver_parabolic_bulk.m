function u = solver_parabolic_bulk(D,f,P,M,K,R,bcond,T,tau,u0)
    
    n = length(D);

    bulknodes = (1:length(M))';
    if strcmp(bcond, 'dir')
        boundarynodes = sum(R,2) == 1;
    else
        boundarynodes = [];
    end
    bulknodes(boundarynodes) = [];

    Mbcond = M(bulknodes, bulknodes);
    Kbcond = K(bulknodes, bulknodes);
    
    for i=1:n
        LHS(:,:,i) = tau*D(i)*Kbcond + Mbcond; %#ok
        perm(:,i) = symamd(LHS(:,:,i)); %#ok
        [L(:,:,i),U(:,:,i)] = lu(LHS(perm(:,i),perm(:,i),i),'vector'); %#ok
    end

    NT = ceil(T/tau);
    u = u0;

    for i=0:NT-1
        for j=1:n
            RHS(:,j) = M*(u(:,j) + tau*f{j}(u, P, i*tau)); %#ok
        end
        RHS = RHS(bulknodes,:);
        for j=1:n
            ubcond(:,j) =  U(:,:,j)\(L(:,:,1)\RHS(perm(:,j),j)); %#ok
            ubcond(perm(:,j),j) = ubcond(:,j); %#ok
        end
        u = zeros(length(M),n);
        u(bulknodes,:) = ubcond;
    end

end