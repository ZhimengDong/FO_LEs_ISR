function [LEs_history, time_axis] = run_ISR_LEs( ...
    dim, state_rhs, augmented_rhs, T_trans, T_sim, dt_orth, ...
    params, nodes, weights)
%RUN_ISR_LES Compute the augmented Lyapunov spectrum with periodic QR steps.

    nodes = nodes(:);
    weights = weights(:);
    N = numel(nodes);
    total_dim = dim * N;

    if numel(weights) ~= N
        error('nodes and weights must have the same length.');
    end

    if isfield(params, 'internal_initial')
        Z0 = params.internal_initial(:);
    elseif isfield(params, 'initial_state')
        Z0 = zeros(total_dim, 1);
    else
        Z0 = 0.1 * rand(total_dim, 1);
    end

    if numel(Z0) ~= total_dim
        error('The augmented initial state must contain dim*N elements.');
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
            @(t, Z) state_rhs(t, Z, params, nodes, weights, N), ...
            [0, T_trans], Z0, ode_options);
        Z_current = Z_transient(end, :)';
    else
        Z_current = Z0;
    end

    n_steps = floor(T_sim / dt_orth);
    Phi_current = eye(total_dim);
    LE_sum = zeros(total_dim, 1);
    LEs_history = zeros(n_steps, total_dim);
    time_axis = (1:n_steps)' * dt_orth;

    for k = 1:n_steps
        Y0 = [Z_current; Phi_current(:)];
        [~, Y_interval] = ode15s( ...
            @(t, Y) augmented_rhs( ...
            t, Y, params, nodes, weights, N, total_dim), ...
            [0, dt_orth], Y0, ode_options);

        Y_end = Y_interval(end, :)';
        Z_current = Y_end(1:total_dim);
        Phi_raw = reshape( ...
            Y_end(total_dim+1:end), total_dim, total_dim);

        [Phi_current, R] = qr(Phi_raw, 0);

        negative_diagonal = diag(R) < 0;
        Phi_current(:, negative_diagonal) = ...
            -Phi_current(:, negative_diagonal);
        R(negative_diagonal, :) = -R(negative_diagonal, :);

        LE_sum = LE_sum + log(abs(diag(R)));
        LEs_history(k, :) = (LE_sum / (k * dt_orth))';
    end
end

