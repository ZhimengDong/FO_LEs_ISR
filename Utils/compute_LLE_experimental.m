function [t_span, log_distance, LLE_estimate, fit_indices, fit_coeff] = ...
    compute_LLE_experimental( ...
    state_rhs, dim, T_trans, T_sim, dt_sample, params, ...
    nodes, weights, fit_time_range)
%COMPUTE_LLE_EXPERIMENTAL Estimate the LLE from two adjacent trajectories.

    nodes = nodes(:);
    weights = weights(:);
    N = numel(nodes);
    total_dim = dim * N;

    if isfield(params, 'internal_initial')
        Z_initial = params.internal_initial(:);
    elseif isfield(params, 'initial_state')
        Z_initial = zeros(total_dim, 1);
    else
        Z_initial = rand(total_dim, 1);
    end

    transient_options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);
    if T_trans > 0
        [~, Z_transient] = ode15s( ...
            @(t, Z) state_rhs(t, Z, params, nodes, weights, N), ...
            [0, T_trans], Z_initial, transient_options);
        Z_reference0 = Z_transient(end, :)';
    else
        Z_reference0 = Z_initial;
    end

    perturbation_size = 1e-10;
    Z_perturbed0 = Z_reference0;
    Z_perturbed0(1) = Z_perturbed0(1) + perturbation_size;

    tracking_options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
    t_span = (0:dt_sample:T_sim)';

    [~, Z_reference] = ode15s( ...
        @(t, Z) state_rhs(t, Z, params, nodes, weights, N), ...
        t_span, Z_reference0, tracking_options);
    [~, Z_perturbed] = ode15s( ...
        @(t, Z) state_rhs(t, Z, params, nodes, weights, N), ...
        t_span, Z_perturbed0, tracking_options);

    distance_squared = zeros(numel(t_span), 1);
    for i_dim = 1:dim
        idx = (i_dim - 1) * N + (1:N);
        x_reference = Z_reference(:, idx) * weights;
        x_perturbed = Z_perturbed(:, idx) * weights;
        distance_squared = distance_squared + ...
            (x_perturbed - x_reference).^2;
    end

    distance = sqrt(distance_squared);
    distance(distance < 1e-16) = 1e-16;
    log_distance = log(distance);

    fit_indices = find( ...
        t_span >= fit_time_range(1) & t_span <= fit_time_range(2));
    fit_coeff = polyfit( ...
        t_span(fit_indices), log_distance(fit_indices), 1);
    LLE_estimate = fit_coeff(1);
end

