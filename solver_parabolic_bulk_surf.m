function [u, v, t, uprime_norm, vprime_norm, u_average, v_average] = solver_parabolic_bulk_surf(dOmega, dGamma, f, g, h, P, M, MS, K, KS, R, T, tau, u0, v0)

n = length(dOmega);
m = length(dGamma);
PGamma = R'*P;

for i=1:n
    LHS_Omega{i} = M+tau*dOmega(i)*K; %#ok
    perm_Omega(:,i) = symamd(LHS_Omega{i}); %#ok
    [L_Omega{i},U_Omega{i}] = lu(LHS_Omega{i}(perm_Omega(:,i),perm_Omega(:,i)),'vector'); %#ok
end

for i=1:m
    LHS_Gamma{i} = MS+tau*dGamma(i)*KS; %#ok
    perm_Gamma(:,i) = symamd(LHS_Gamma{i}); %#ok
    [L_Gamma{i},U_Gamma{i}] = lu(LHS_Gamma{i}(perm_Gamma(:,i),perm_Gamma(:,i)),'vector'); %#ok
end

NT = ceil(T/tau);
u =  u0;
v =  v0;
if nargout >= 4
     % spatial L2 norm of time derivative of fist component of u and v
     uprime_norm = zeros(1,NT);
     vprime_norm = zeros(1,NT);
end
if nargout >= 6
     % spatial average of first component of u and v
     u_average = zeros(1,NT);
     v_average = zeros(1,NT);
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
       RHS_Omega(:,j)= M*(u(:,j) + tau*f{j}(u,P,i*tau)) + dOmega(j)*tau*R*MS*h{j}(R'*u,v, PGamma, i*tau); %#ok
   end
   for j=1:m
       RHS_Gamma(:,j)= MS*(v(:,j) + tau*g{j}(R'*u,v,PGamma,i*tau)); %#ok
   end

   for j=1:n
       unew(:,j) =  U_Omega{j}\(L_Omega{j}\RHS_Omega(perm_Omega(:,j),j)); %#ok
       unew(perm_Omega(:,j),j) = unew(:,j); %#ok
   end
   for j=1:m
       vnew(:,j) =  U_Gamma{j}\(L_Gamma{j}\RHS_Gamma(perm_Gamma(:,j),j)); %#ok
       vnew(perm_Gamma(:,j),j) = vnew(:,j); %#ok
   end

   if nargout >= 4
       incr_u = unew(:,1)-u(:,1);
       incr_v = vnew(:,1)-v(:,1);
       uprime_norm(i+1) = incr_u'*M*incr_u;
       vprime_norm(i+1) = incr_v'*MS*incr_v;
   end
   if nargout >= 6
       u_average(i+1) = sum(M*unew(:,1));
       v_average(i+1) = sum(MS*vnew(:,1));
   end

   u = unew;
   v = vnew;

end

close(progress_handle);

if nargout >= 3
    t = linspace(tau, NT*tau, NT);
end
if nargout >= 4
    uprime_norm = sqrt(uprime_norm);
    vprime_norm = sqrt(vprime_norm);
end
if nargout >= 6
    u_average = u_average/sum(sum(M));
    v_average = v_average/sum(sum(MS));
end

end