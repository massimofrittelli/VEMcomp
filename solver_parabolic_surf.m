function [v, t, vprime_norm, v_average] = solver_parabolic_surf(D,g,P,MS,KS,R,T,tau,v0)
    
    n = length(D);

    for i=1:n
        LHS{i} = tau*D(i)*KS + MS; %#ok
        perm(:,i) = symamd(LHS{i}); %#ok
        [L{i},U{i}] = lu(LHS{i}(perm(:,i),perm(:,i)),'vector'); %#ok
    end
    PS = R'*P;
    NT = ceil(T/tau);

    v = v0;
    v_new = zeros(size(v0));
    if nargout >= 3
        % spatial L2 norm of time derivative of fist component of v
        vprime_norm = zeros(1,NT);
    end
    if nargout >= 4
        % spatial average of first component of v
        v_average = zeros(1,NT);
    end
    for i=0:NT-1
        for j=1:n
            RHS(:,j) = MS*(v(:,j) + tau*g{j}(v, PS, i*tau)); %#ok
            v_new(perm(:,j),j) =  U{j}\(L{j}\RHS(perm(:,j),j));
        end
        if nargout >= 3
            incr = v_new(:,1)-v(:,1);
            vprime_norm(i+1) = incr'*MS*incr;
        end
        if nargout >= 4
            v_average(i+1) = sum(MS*v_new(:,1));
        end
        v = v_new;
    end

    if nargout >= 2
        t = linspace(tau, NT*tau, NT);
    end
    if nargout >= 3
        vprime_norm = sqrt(vprime_norm);
    end
    if nargout >= 4
        v_average = v_average/sum(sum(MS));
    end

end