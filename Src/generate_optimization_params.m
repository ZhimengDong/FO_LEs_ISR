function [omega_opt, eta_opt] = generate_optimization_params( ...
    alpha, N, dt, T_memory)
%GENERATE_OPTIMIZATION_PARAMS Fit the fractional kernel in the time domain.
%   This optional alternative requires Optimization Toolbox.

    target = @(t) t.^(alpha - 1) ./ gamma(alpha);
    t_samples = logspace(log10(dt), log10(T_memory), 200)';
    y_target = target(t_samples);

    omega0 = logspace(log10(1 / T_memory), log10(1 / dt), N);
    eta0 = ones(1, N) * (mean(y_target) / N);
    x0 = [eta0, omega0];

    lower_bound = ones(1, 2 * N) * 1e-9;
    upper_bound = inf(1, 2 * N);

    options = optimoptions( ...
        'lsqcurvefit', ...
        'Display', 'off', ...
        'MaxFunctionEvaluations', 5000, ...
        'StepTolerance', 1e-8, ...
        'FunctionTolerance', 1e-8, ...
        'Algorithm', 'trust-region-reflective');

    model = @(x, t) multi_exponential_kernel(x, t, N);
    x_opt = lsqcurvefit( ...
        model, x0, t_samples, y_target, ...
        lower_bound, upper_bound, options);

    eta_opt = x_opt(1:N)';
    omega_opt = x_opt(N+1:end)';
end

function y = multi_exponential_kernel(x, t, N)
    eta = x(1:N);
    omega = x(N+1:end);
    y = exp(-t * omega) * eta';
end

