function [LEs_history, Time_axis] = run_ISR_LEs(dim, sys_deriv, augmented_sys, ...
    T_trans, T_sim, dt_orth, params, nodes, weights)
    % RUN_ISR_LES Computes the complete Lyapunov spectrum of the augmented proxy system.
    %
    % Inputs:
    %   dim  : Dimension of the original physical system (e.g., 3 for Lorenz)
    %   sys_deriv  : Function handle for the augmented state equations
    %   augmented_sys  : Function handle for the combined state + variational equations
    %   T_trans  : Time to discard transients
    %   T_sim  : Total simulation time for LEs computation
    %   dt_orth  : Time interval for periodic Gram-Schmidt reorthonormalization
    %   params  : Struct containing system parameters
    %   nodes, weights  : Frequency nodes and weights from the approximation kernel
    %
    % Outputs:
    %   LEs_history  : Matrix of shape [n_steps, total_dim] tracking LEs convergence
    %   Time_axis  : Time vector corresponding to the steps

    dim_per_var = length(nodes);
    total_dim = dim * dim_per_var;  % Augmented state space dimension M = n * Q
    
    % 1. Eliminate transients to ensure the trajectory is on the attractor
    Z0 = 0.1 * rand(total_dim, 1);  % Random initial value
    opt_ode = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);  % Tolerence error
    
    [~, Z_traj] = ode15s(@(t,Z) sys_deriv(t, Z, params, nodes, weights, dim_per_var), ...
                         [0, T_trans], Z0, opt_ode);
    Z_start = Z_traj(end, :)';
    
    % 2. Lyapunov spectrum computation loop
    n_steps = floor(T_sim / dt_orth);
    
    Z_curr = Z_start;
    Q_curr = eye(total_dim);  % Initialize orthonormal basis
    LE_sum = zeros(total_dim, 1);
    LEs_history = zeros(n_steps, total_dim);
    Time_axis = (1:n_steps) * dt_orth;
    
    for k = 1:n_steps
        % Flatten the state and variational matrix for the ODE solver
        Y0 = [Z_curr; Q_curr(:)];
        
        % Integrate over one orthonormalization interval
        [~, Y_res] = ode15s(@(t,Y) augmented_sys(t, Y, params, nodes, weights, dim_per_var, total_dim), ...
                            [0, dt_orth], Y0, opt_ode);
        Y_end = Y_res(end, :)';
        
        % Extract evolved state and deformed variational matrix
        Z_next = Y_end(1:total_dim);
        Q_raw = reshape(Y_end(total_dim+1:end), total_dim, total_dim);
        
        % Perform discrete reorthonormalization (via QR decomposition)
        [Q_next, R] = qr(Q_raw, 0);
        
        % Ensure positive diagonal elements in R for continuous basis evolution
        idx_neg = find(diag(R) < 0);
        Q_next(:, idx_neg) = -Q_next(:, idx_neg);
        R(idx_neg, :) = -R(idx_neg, :);
        
        % Accumulate local expansion rates
        LE_sum = LE_sum + log(abs(diag(R)));
        LEs_history(k, :) = LE_sum' / (k * dt_orth);
        
        % Reset for the next iteration vaildly
        Z_curr = Z_next;
        Q_curr = Q_next;
    end
end