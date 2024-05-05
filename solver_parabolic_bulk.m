function [u, t, uprime_norm, u_average] = solver_parabolic_bulk(D,f,P,M,K,R,bcond,T,tau,u0)
    
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
        LHS{i} = tau*D(i)*Kbcond + Mbcond; %#ok
        perm(:,i) = symamd(LHS{i}); %#ok
        [L{i},U{i}] = lu(LHS{i}(perm(:,i),perm(:,i)),'vector'); %#ok
    end

    NT = ceil(T/tau);
    u = u0;
    if nargout >= 3
        % spatial L2 norm of time derivative of fist component of u
        uprime_norm = zeros(1,NT);
    end
    if nargout >= 4
        % spatial average of first component of u
        u_average = zeros(1,NT);
    end
    progress_handle = waitbar(0);
    axes_handle = findobj(progress_handle, 'type','axes');
    title_handle = get(axes_handle,'title');
    set(title_handle, 'FontSize', 18);
    waitbar(0, progress_handle, 'Timestepping in progress: 0 %')

    percent_prev = 0;
    for i=0:NT-1
        percent_new = round(i*100/NT);
        if percent_new > percent_prev
            waitbar(i/NT,progress_handle, sprintf('Timestepping in progress: %d%%', percent_new));
            percent_prev = percent_new;
        end
        for j=1:n
            RHS(:,j) = M*(u(:,j) + tau*f{j}(u, P, i*tau)); %#ok
        end
        RHS = RHS(bulknodes,:);
        for j=1:n
            ubcond(:,j) =  U{j}\(L{j}\RHS(perm(:,j),j)); %#ok
            ubcond(perm(:,j),j) = ubcond(:,j); %#ok
        end
        u_new = zeros(length(M),n);
        u_new(bulknodes,:) = ubcond;
        if nargout >= 3
            incr = u_new(:,1)-u(:,1);
            uprime_norm(i+1) = incr'*M*incr;
        end
        if nargout >= 4
            u_average(i+1) = sum(M*u_new(:,1));
        end
        u = u_new;
    end

    close(progress_handle);

    if nargout >= 2
        t = linspace(tau, NT*tau, NT);
    end
    if nargout >= 3
        uprime_norm = sqrt(uprime_norm);
    end
    if nargout >= 4
        u_average = u_average/sum(sum(M));
    end
end