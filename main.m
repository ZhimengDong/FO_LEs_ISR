function main(recompute)
%% Fractional Lorenz Lyapunov exponents using the proposed framework
% Example corresponding to Fig. 4 of the paper.

clc;
close all;

addpath(fullfile(pwd, 'Src'), fullfile(pwd, 'Utils'));

if nargin < 1
    recompute = false;
end

paper_data_path = fullfile('Data', 'Lorenz_ISR_Results_alpha_0.995.mat');
delta = 1e-1;

if recompute
    %% System and numerical settings
    params.sigma = 10;
    params.rho = 28;
    params.beta = 8/3;
    params.alpha = 0.995;
    params.ode_rel_tol = 1e-10;
    params.ode_abs_tol = 1e-10;

    dim = 3;
    T_trans = 50;
    T_sim = 200;
    dt_orth = 0.5;
    dt_sample = 0.05;
    fit_time_range = [0, 10];

    N = 10;
    w_low = 1e-4;
    w_high = 1e4;

    fprintf('Generating the Oustaloup approximation...\n');
    [nodes, weights] = generate_oustaloup_params( ...
        params.alpha, N - 1, w_low, w_high);

    params.internal_initial = 0.1 * rand(dim * N, 1);

    fprintf('Computing the physical time series...\n');
    [t_seq, X_seq, Y_seq, Z_seq] = generate_ISR_timeseries( ...
        @lorenz_isr_state, dim, T_trans, T_sim, dt_sample, ...
        params, nodes, weights);

    fprintf('Computing the augmented Lyapunov spectrum...\n');
    tic;
    [LEs_history, Time] = run_ISR_LEs( ...
        dim, @lorenz_isr_state, @lorenz_isr_augmented, ...
        T_trans, T_sim, dt_orth, params, nodes, weights);
    runtime_seconds = toc;

    fprintf('Computing the two-trajectory LLE estimate...\n');
    [t_exp, log_dist, LLE_value, idx_fit, P_fit] = ...
        compute_LLE_experimental( ...
        @lorenz_isr_state, dim, T_trans, T_sim, 0.02, ...
        params, nodes, weights, fit_time_range);

    Lorenz_ISR_Results = struct();
    Lorenz_ISR_Results.t_seq = t_seq;
    Lorenz_ISR_Results.X_seq = X_seq;
    Lorenz_ISR_Results.Y_seq = Y_seq;
    Lorenz_ISR_Results.Z_seq = Z_seq;
    Lorenz_ISR_Results.Time = Time;
    Lorenz_ISR_Results.LEs_history = LEs_history;
    Lorenz_ISR_Results.t_exp = t_exp;
    Lorenz_ISR_Results.log_dist = log_dist;
    Lorenz_ISR_Results.LLE_value = LLE_value;
    Lorenz_ISR_Results.idx_fit = idx_fit;
    Lorenz_ISR_Results.P_fit = P_fit;
    Lorenz_ISR_Results.runtime_seconds = runtime_seconds;
    Lorenz_ISR_Results.settings = struct( ...
        'alpha', params.alpha, 'N', N, 'w_low', w_low, ...
        'w_high', w_high, 'T_trans', T_trans, 'T_sim', T_sim, ...
        'dt_orth', dt_orth, 'delta', delta);

    result_path = fullfile( ...
        'Data', 'Lorenz_ISR_Results_alpha_0.995_recomputed.mat');
    save(result_path, 'Lorenz_ISR_Results', '-v7.3');
    output_path = 'Lorenz_ISR_Results_recomputed.png';
else
    if ~isfile(paper_data_path)
        error('Paper data not found: %s', paper_data_path);
    end
    result_path = paper_data_path;
    output_path = 'Lorenz_ISR_Results.png';
end

plot_Lorenz_ISR_results(result_path, output_path, delta);
end

%% Local Lorenz system functions
function dZ = lorenz_isr_state(~, Z, p, nodes, weights, N)
    Zx = Z(1:N);
    Zy = Z(N+1:2*N);
    Zz = Z(2*N+1:3*N);

    x0 = zeros(3, 1);
    if isfield(p, 'initial_state')
        x0 = p.initial_state(:);
    end

    x = x0(1) + weights(:)' * Zx;
    y = x0(2) + weights(:)' * Zy;
    z = x0(3) + weights(:)' * Zz;

    fx = p.sigma * (y - x);
    fy = p.rho * x - x * z - y;
    fz = x * y - p.beta * z;

    dZ = [-nodes(:) .* Zx + fx; ...
          -nodes(:) .* Zy + fy; ...
          -nodes(:) .* Zz + fz];
end

function dY = lorenz_isr_augmented(~, Y, p, nodes, weights, N, total_dim)
    Z = Y(1:total_dim);
    Phi = reshape(Y(total_dim+1:end), total_dim, total_dim);

    Zx = Z(1:N);
    Zy = Z(N+1:2*N);
    Zz = Z(2*N+1:3*N);

    x0 = zeros(3, 1);
    if isfield(p, 'initial_state')
        x0 = p.initial_state(:);
    end

    x = x0(1) + weights(:)' * Zx;
    y = x0(2) + weights(:)' * Zy;
    z = x0(3) + weights(:)' * Zz;

    fx = p.sigma * (y - x);
    fy = p.rho * x - x * z - y;
    fz = x * y - p.beta * z;
    dZ = [-nodes(:) .* Zx + fx; ...
          -nodes(:) .* Zy + fy; ...
          -nodes(:) .* Zz + fz];

    C = repmat(weights(:)', N, 1);
    J_xx = diag(-nodes) - p.sigma * C;
    J_xy = p.sigma * C;
    J_xz = zeros(N);

    J_yx = (p.rho - z) * C;
    J_yy = diag(-nodes) - C;
    J_yz = -x * C;

    J_zx = y * C;
    J_zy = x * C;
    J_zz = diag(-nodes) - p.beta * C;

    J_aug = [J_xx, J_xy, J_xz; ...
             J_yx, J_yy, J_yz; ...
             J_zx, J_zy, J_zz];

    dPhi = J_aug * Phi;
    dY = [dZ; dPhi(:)];
end
