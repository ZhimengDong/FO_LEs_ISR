function [t_out, X, Y, Z, W] = generate_ISR_timeseries( ...
    state_rhs, dim, T_trans, T_sim, dt_sample, params, nodes, weights)
%GENERATE_ISR_TIMESERIES Generate uniformly sampled physical trajectories.

    nodes = nodes(:);
    weights = weights(:);
    N = numel(nodes);
    total_dim = dim * N;

    if isfield(params, 'internal_initial')
        Z0 = params.internal_initial(:);
    elseif isfield(params, 'initial_state')
        Z0 = zeros(total_dim, 1);
    else
        Z0 = 0.1 * rand(total_dim, 1);
    end

    rel_tol = 1e-10;
    abs_tol = 1e-10;
    if isfield(params, 'ode_rel_tol')
        rel_tol = params.ode_rel_tol;
    end
    if isfield(params, 'ode_abs_tol')
        abs_tol = params.ode_abs_tol;
    end
    ode_options = odeset('RelTol', rel_tol, 'AbsTol', abs_tol);

    if T_trans > 0
        [~, Z_transient] = ode15s( ...
            @(t, Zs) state_rhs(t, Zs, params, nodes, weights, N), ...
            [0, T_trans], Z0, ode_options);
        Z_start = Z_transient(end, :)';
    else
        Z_start = Z0;
    end

    sample_times = (0:dt_sample:T_sim)';
    [t_out, Z_trajectory] = ode15s( ...
        @(t, Zs) state_rhs(t, Zs, params, nodes, weights, N), ...
        sample_times, Z_start, ode_options);

    x0 = zeros(dim, 1);
    if isfield(params, 'initial_state')
        x0 = params.initial_state(:);
    end

    physical = zeros(numel(t_out), dim);
    for i_dim = 1:dim
        idx = (i_dim - 1) * N + (1:N);
        physical(:, i_dim) = ...
            x0(i_dim) + Z_trajectory(:, idx) * weights;
    end

    X = physical(:, 1);
    Y = [];
    Z = [];
    W = [];
    if dim >= 2
        Y = physical(:, 2);
    end
    if dim >= 3
        Z = physical(:, 3);
    end
    if dim >= 4
        W = physical(:, 4);
    end
end
