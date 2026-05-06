function [t_out, X, Y, Z] = generate_ISR_timeseries(sys_deriv, dim, T_trans, T_sim, dt_sample, params, nodes, weights)
    % GENERATE_ISR_TIMESERIES Extracts fixed-step physical time series from the ISR framework.
    % Essential for phase space reconstruction and the 0-1 test for chaos.
    
    Q = length(nodes);
    total_dim = dim * Q;
    
    % Initial conditions
    Z0 = 0.1 * rand(total_dim, 1);
    opt_ode = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);
    
    % 1. Integrate out the transient behavior
    [~, Z_trans] = ode15s(@(t,Z) sys_deriv(t, Z, params, nodes, weights, Q), ...
                         [0, T_trans], Z0, opt_ode);
    Z_start = Z_trans(end, :)';
    
    % 2. Formal integration with fixed sampling steps
    t_span = 0 : dt_sample : T_sim;
    [t_out, Z_traj] = ode15s(@(t,Z) sys_deriv(t, Z, params, nodes, weights, Q), ...
                             t_span, Z_start, opt_ode);
    
    % 3. Project the augmented high-dimensional state back to physical phase space
    if dim == 3
        Zx = Z_traj(:, 1:Q);
        Zy = Z_traj(:, Q+1:2*Q);
        Zz = Z_traj(:, 2*Q+1:3*Q);
        
        X = Zx * weights(:);
        Y = Zy * weights(:);
        Z = Zz * weights(:);
    else
        Zx = Z_traj(:, 1:Q);
        Zy = Z_traj(:, Q+1:2*Q);

        X = Zx * weights(:);
        Y = Zy * weights(:);
        Z = [];
    end
end