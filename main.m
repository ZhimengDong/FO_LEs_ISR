%% ========================================================================
% Main Script: Computation of Fractional LEs using ISR Framework
% Corresponding to the paper:
% [Mitigating algorithmic false chaos in fractional-order systems via infinite state representation]
%
% This script computes the complete Lyapunov Exponent (LE) spectrum, 
% evaluates the 0-1 test, and performs the two-trajectory experimental 
% LLE validation for the fractional-order Lorenz system.
% ========================================================================

clc; clear; close all;
addpath('./Src', './Utils');  % Ensure subfolders are in the path

%% 1. System and Simulation Parameters
sys_p.sigma = 10;
sys_p.rho = 28;
sys_p.beta = 8/3;
sys_p.alpha = 0.97;  % Fractional order

dim = 3;  % Dimension of the original physical system
T_trans = 50;  % Transient time to discard
T_sim = 200;  % Total simulation time for LEs evolution
dt_orth = 0.5;  % Time interval for Gram-Schmidt reorthonormalization
dt_sample = 0.05;  % Sampling step for time series generation
fit_time_range = [0, 10];  % Time window for experimental LLE linear fitting

%% 2. Generate Frequency Nodes and Weights (ISR Approximation)
Q = 10;  % Number of frequency nodes per physical dimension

% Frequency band limits for Oustaloup / Time limits for Optimization
w_low  = 1e-4;  
w_high = 1e4;   
dt_fit = 1/w_high;   
T_mem  = 1/w_low;

fprintf('Generating Approximation Kernel (Method A: Oustaloup)...\n');
[nodes, weights] = generate_oustaloup_params(sys_p.alpha, Q-1, w_low, w_high);

% Method B: Time-domain Optimization (Section 5.3)
% fprintf('Generating Approximation Kernel (Method B: Optimization)...\n');
% [nodes, weights] = generate_optimization_params(sys_p.alpha, Q, dt_fit, T_mem);

%% 3. Execute Computations
fprintf('Computing continuous time series...\n');
[t_seq, X_seq, Y_seq, Z_seq] = generate_ISR_timeseries(@sys_deriv_Lorenz, dim, ...
    T_trans, T_sim, dt_sample, sys_p, nodes, weights);

tic;
fprintf('\nComputing the complete augmented Lyapunov spectrum...\n');
[LEs_history, Time] = run_ISR_LEs(dim, @sys_deriv_Lorenz, @augmented_sys_Lorenz, ...
    T_trans, T_sim, dt_orth, sys_p, nodes, weights);
toc;

fprintf('\nCross-validating via Experimental LLE (Two-trajectory method)...\n');
[t_exp, log_dist, LLE_value, idx_fit, P_fit] = compute_LLE_experimental(...
    @sys_deriv_Lorenz, dim, T_trans, T_sim, 0.02, sys_p, nodes, weights, fit_time_range);

fprintf('\nEvaluation complete. LLE extracted: %.4f\n', LLE_value);

%% 4. Save Results and Plot
% Ensure the directory exists before saving
if ~exist('./Data', 'dir'), mkdir('./Data'); end

Lorenz_ISR_Results.X_seq = X_seq;
Lorenz_ISR_Results.Time = Time;
Lorenz_ISR_Results.LEs_history = LEs_history;
Lorenz_ISR_Results.t_exp = t_exp;
Lorenz_ISR_Results.log_dist = log_dist;
Lorenz_ISR_Results.LLE_value = LLE_value;
Lorenz_ISR_Results.idx_fit = idx_fit;
Lorenz_ISR_Results.P_fit = P_fit;

save_path = fullfile('./Data', 'Lorenz_ISR_Results_alpha_0.97.mat');
save(save_path, 'Lorenz_ISR_Results', '-v7.3');
fprintf('Results saved to %s\n', save_path);

% visualization
fprintf('\nGenerating figures...\n');
run('plot_Lorenz_ISR_results.m');


%% ========================================================================
% LOCAL FUNCTIONS: System Definitions
% ========================================================================

% State equations for the augmented fractional Lorenz system
function dZ = sys_deriv_Lorenz(~, Z, p, w, c, dim)
    % Unpack augmented states
    Zx = Z(1:dim); Zy = Z(dim+1:2*dim); Zz = Z(2*dim+1:3*dim);
    
    % Project memory variables back to physical variables: X(t) = sum(c_i * Zx_i)
    X = c' * Zx; 
    Y = c' * Zy; 
    Z_val = c' * Zz;  % New name
    
    % Original Lorenz nonlinear vector field
    fx = p.sigma * (Y - X);
    fy = p.rho * X - X * Z_val - Y;
    fz = X * Y - p.beta * Z_val;
    
    % Distributed state derivatives: dZ_i = -w_i * Z_i + f_original
    dZx = -w .* Zx + fx;
    dZy = -w .* Zy + fy;
    dZz = -w .* Zz + fz;
    
    dZ = [dZx; dZy; dZz];
end

% Augmented variational equations (State + Jacobian evolution)
function dY = augmented_sys_Lorenz(~, Y, p, w, c, dim, total_dim)
    Z = Y(1:total_dim);
    Q_mat = reshape(Y(total_dim+1:end), total_dim, total_dim);
    
    Zx = Z(1:dim); Zy = Z(dim+1:2*dim); Zz = Z(2*dim+1:3*dim);
    X = c' * Zx; Y = c' * Zy; Z_val = c' * Zz;
    
    % 1. State derivatives
    fx = p.sigma * (Y - X);
    fy = p.rho * X - X * Z_val - Y;
    fz = X * Y - p.beta * Z_val;
    dZ = [-w .* Zx + fx; -w .* Zy + fy; -w .* Zz + fz];
    
    % 2. Construct Augmented Jacobian Matrix J (total_dim x total_dim)
    % Precompute weight matrix (contribution of all modes to nonlinear inputs is c_i)
    C_mat = repmat(c', dim, 1); 
    
    % Block 1: d(dZx)/dZ
    J_xx = diag(-w) - p.sigma * C_mat;
    J_xy = p.sigma * C_mat;
    J_xz = zeros(dim);
    
    % Block 2: d(dZy)/dZ
    J_yx = (p.rho - Z_val) * C_mat;
    J_yy = diag(-w) - C_mat;
    J_yz = -X * C_mat;
    
    % Block 3: d(dZz)/dZ
    J_zx = Y * C_mat;
    J_zy = X * C_mat;
    J_zz = diag(-w) - p.beta * C_mat;
    
    Jac = [J_xx, J_xy, J_xz;
           J_yx, J_yy, J_yz;
           J_zx, J_zy, J_zz];
           
    % 3. Variational matrix evolution: dPhi = J * Phi
    dQ = Jac * Q_mat;
    
    dY = [dZ; dQ(:)];
end
