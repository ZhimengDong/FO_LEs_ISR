function [t_span, log_dist, LLE_est, idx_fit, P_fit] = compute_LLE_experimental(...
    sys_deriv, dim, T_trans, T_sim, dt_sample, params, nodes, weights, fit_time_range)
    % COMPUTE_LLE_EXPERIMENTAL Estimates the Largest Lyapunov Exponent (LLE)
    % using the classical two-trajectory tracking method.
    % Acts as an independent cross-validation.
    
    Q = length(nodes);
    total_dim = dim * Q;
    
    % 1. Transient phase to ensure initial points strictly reside on the attractor
    Z_init = rand(total_dim, 1);
    opt_trans = odeset('RelTol', 1e-8, 'AbsTol', 1e-8); 
    [~, Z_trans] = ode15s(@(t,Z) sys_deriv(t, Z, params, nodes, weights, Q), ...
                         [0, T_trans], Z_init, opt_trans);
    Z0_ref = Z_trans(end, :)';
    
    % 2. Apply infinitesimal perturbation to the initial state
    delta0 = 1e-10; 
    Z0_pert = Z0_ref;
    % Perturb the first frequency mode of the first physical variable
    Z0_pert(1) = Z0_pert(1) + delta0; 
    
    % 3. Synchronously integrate the reference and perturbed trajectories
    % High solver precision is mandatory to prevent numerical noise from swamping the perturbation
    opt_track = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
    t_span = (0:dt_sample:T_sim)';
    
    [~, Z_ref_traj]  = ode15s(@(t,Z) sys_deriv(t, Z, params, nodes, weights, Q), t_span, Z0_ref, opt_track);
    [~, Z_pert_traj] = ode15s(@(t,Z) sys_deriv(t, Z, params, nodes, weights, Q), t_span, Z0_pert, opt_track);
    
    % 4. Compute physical distance between the two trajectories
    if dim == 3
        X_ref = Z_ref_traj(:, 1:Q) * weights(:); Y_ref = Z_ref_traj(:, Q+1:2*Q) * weights(:); Z_ref = Z_ref_traj(:, 2*Q+1:3*Q) * weights(:);
        X_pert = Z_pert_traj(:, 1:Q) * weights(:); Y_pert = Z_pert_traj(:, Q+1:2*Q) * weights(:); Z_pert = Z_pert_traj(:, 2*Q+1:3*Q) * weights(:);
        
        dist = sqrt((X_pert - X_ref).^2 + (Y_pert - Y_ref).^2 + (Z_pert - Z_ref).^2);
    else
        X_ref = Z_ref_traj(:, 1:Q) * weights(:); Y_ref = Z_ref_traj(:, Q+1:2*Q) * weights(:);
        X_pert = Z_pert_traj(:, 1:Q) * weights(:); Y_pert = Z_pert_traj(:, Q+1:2*Q) * weights(:);
        
        dist = sqrt((X_pert - X_ref).^2 + (Y_pert - Y_ref).^2);
    end
    
    % Logarithmic transform.
    dist(dist < 1e-16) = 1e-16; 
    log_dist = log(dist);
    
    % 5. Linear regression in the designated exponential growth phase to extract LLE
    idx_fit = find(t_span >= fit_time_range(1) & t_span <= fit_time_range(2));
    P_fit = polyfit(t_span(idx_fit), log_dist(idx_fit), 1);
    LLE_est = P_fit(1);
end