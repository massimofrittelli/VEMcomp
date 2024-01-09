function [L2_relative_err] = compute_error(C,MS,u,u_exact,v,v_exact)
    
    if size(u) ~= size(u_exact)
        error('u and u_exact must have the same size!')
    end

    if size(v) ~= size(v_exact)
        error('v and v_exact must have the same size!')
    end

    pw_err_bulk = u - u_exact;
    pw_err_surf = v - v_exact;

    square_norm_sol_bulk = zeros(1, size(u,2));
    square_err_bulk = zeros(1, size(u,2));
    square_norm_sol_surf = zeros(1, size(v,2));
    square_err_surf = zeros(1, size(v,2));

    for i=1:size(u,2)
        square_norm_sol_bulk(i) = u_exact(:,i)'*C*u_exact(:,i);
        square_err_bulk(i) = pw_err_bulk(:,i)'*C*pw_err_bulk(:,i);
    end
    
    for i=1:size(v,2)
        square_norm_sol_surf(i) = v_exact(:,i)'*MS*v_exact(:,i);
        square_err_surf(i) = pw_err_surf(:,i)'*MS*pw_err_surf(:,i);
    end

    L2_relative_err = sqrt((sum(square_err_bulk) + sum(square_err_surf)) ...
        /(sum(square_norm_sol_bulk) + sum(square_norm_sol_surf)));

end