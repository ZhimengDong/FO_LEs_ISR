function [omega_opt, w_opt] = generate_optimization_params(alpha, Q, dt, T_mem)
    % GENERATE_OPTIMIZATION_PARAMS Generates nodes and weights via 
    % time-domain nonlinear least-squares optimization. 
    % (Method B in the associated manuscript, used for cross-validation).
    %
    % Objective: Fit the fractional kernel K(t) = sum(w_k * exp(-omega_k * t))
    
    target_func = @(t) (t.^(alpha-1)) ./ gamma(alpha);
    
    % Distribute sample points logarithmically to balance short-term (high freq)
    % and long-term (low freq) memory effects
    t_samples = logspace(log10(dt), log10(T_mem), 200)';
    y_target = target_func(t_samples);
    
    % Initial Guess: Crucial to avoid local minima. 
    % Use logarithmic uniform distribution for initial frequency nodes.
    omega0 = logspace(log10(1/T_mem), log10(1/dt), Q);
    w0 = ones(1, Q) * (mean(y_target)/Q); 
    x0 = [w0, omega0];
    
    % Constraints: Weights and frequencies must be strictly positive
    lb = zeros(1, 2*Q) + 1e-9; 
    ub = inf(1, 2*Q);
    
    % Optimization configurations
    options = optimoptions('lsqcurvefit', 'Display', 'off', ...
        'MaxFunctionEvaluations', 5000, 'StepTolerance', 1e-8, ...
        'FunctionTolerance', 1e-8, 'Algorithm', 'trust-region-reflective');
        
    % Define the multi-exponential fit model
    model_fun = @(x, t) fit_kernel_model(x, t, Q);
    
    % Execute nonlinear least-squares curve fitting
    x_opt = lsqcurvefit(model_fun, x0, t_samples, y_target, lb, ub, options);
    
    w_opt = x_opt(1:Q)';
    omega_opt = x_opt(Q+1:end)';
end

function y = fit_kernel_model(x, t, Q)
    w = x(1:Q);
    omega = x(Q+1:end);
    y = zeros(size(t));
    for k = 1:Q
        y = y + w(k) * exp(-omega(k) * t);
    end
end