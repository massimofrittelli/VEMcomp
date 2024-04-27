function v = solver_parabolic_surf(D,g,P,MS,KS,R,T,tau,v0)
    
    n = length(D);

    for i=1:n
        LHS{i} = tau*D(i)*KS + MS; %#ok
        perm(:,i) = symamd(LHS{i}); %#ok
        [L{i},U{i}] = lu(LHS{i}(perm(:,i),perm(:,i)),'vector'); %#ok
    end
    PS = R'*P;
    NT = ceil(T/tau);

    v = v0;
    for i=0:NT-1
        for j=1:n
            RHS(:,j) = MS*(v(:,j) + tau*g{j}(v, PS, i*tau)); %#ok
            v(perm(:,j),j) =  U{j}\(L{j}\RHS(perm(:,j),j));
        end
    end

end