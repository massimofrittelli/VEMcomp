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
        LHS{i} = tau*D(i)*Kbcond + Mbcond; %#ok
        perm(:,i) = symamd(LHS{i}); %#ok
        [L{i},U{i}] = lu(LHS{i}(perm(:,i),perm(:,i)),'vector'); %#ok
    end

    NT = ceil(T/tau);
    u = u0;
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
        u = zeros(length(M),n);
        u(bulknodes,:) = ubcond;
    end

    close(progress_handle);

end