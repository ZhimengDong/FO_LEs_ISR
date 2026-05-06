%% ========================================================================
% Plotting Script: Visualization of ISR LEs and Cross-Validations
%
% This script loads the pre-computed data and performs:
% (a) Extraction of physical LEs using the Degeneracy-Rejection algorithm.
% (b) Plotting of the LE convergence history.
% (c) Plotting of the Experimental LLE (two-trajectory tracking).
% (d) Plotting of the 0-1 test for chaos (p-q trajectory).
% ========================================================================


%% 1. Load Pre-computed Data
data_dir = './Data';
file_name = 'Lorenz_ISR_Results_alpha_0.97.mat';  % Change alpha value as needed
load_path = fullfile(data_dir, file_name);

if ~exist(load_path, 'file')
    error('Data file not found!');
end
load(load_path);

% Unpack data structures
X_seq = Lorenz_ISR_Results.X_seq;
Time = Lorenz_ISR_Results.Time;
LEs_history = Lorenz_ISR_Results.LEs_history;
LEs = sort(LEs_history(end, :)', 'descend'); % Sorted full spectrum
t_exp = Lorenz_ISR_Results.t_exp;
log_dist = Lorenz_ISR_Results.log_dist;
LLE_value = Lorenz_ISR_Results.LLE_value;
idx_fit = Lorenz_ISR_Results.idx_fit;
P_fit = Lorenz_ISR_Results.P_fit;

%% 2. Degeneracy-Rejection Algorithm (Extracting Physical LEs)
fprintf('Executing Degeneracy-Rejection Algorithm...\n');
n_sys = 3;  % Physical dimension of the Lorenz system
N = length(LEs);
Q = (N - n_sys) / n_sys;  % Number of spurious clusters

% Combinatorial search for the optimal physical LEs extraction
combos = nchoosek(1:N, n_sys);  
min_cost = inf;  
best_phys_idx = [];  
best_spurious_mat = [];  

for i = 1:size(combos, 1)
    phys_idx = combos(i, :);
    
    % Generate mask to separate the assumed physical LEs
    mask = true(N, 1);
    mask(phys_idx) = false;
    spurious_roots = LEs(mask);
    
    % Reshape the remaining (N - n_sys) roots into clusters.
    % Each column represents an n_sys-fold degenerate spurious cluster.
    spurious_mat = reshape(spurious_roots, n_sys, Q);
    
    % Calculate the objective cost function (Relative dispersion)
    cluster_spans = max(spurious_mat) - min(spurious_mat);
    cluster_max_abs = max(abs(spurious_mat));
    
    cost = sum(cluster_spans ./ max(cluster_max_abs, 1e-1));  % Eq. 22
    
    if cost < min_cost
        min_cost = cost;
        best_phys_idx = phys_idx;
        best_spurious_mat = spurious_mat;
    end
end

physical_LEs = LEs(best_phys_idx);  
fprintf('Extraction complete. Physical LEs: [%.4f, %.4f, %.4f]\n', physical_LEs);

%% 3. Execute 0-1 Test for Chaos
fprintf('Computing 0-1 Test for chaos...\n');
[K_val, p_c, q_c] = compute_01_test(X_seq);
fprintf('0-1 Test Asymptotic K-median value: %.4f\n', K_val);

%% 4. Publication-Quality Plotting
figure('Position', [50, 50, 900, 750], 'Color', 'w');

set(0, 'DefaultAxesFontName', 'Times New Roman', 'DefaultAxesFontSize', 12);
set(0, 'DefaultTextFontName', 'Times New Roman', 'DefaultTextFontSize', 12);

color_blue   = [0.00, 0.45, 0.74];
color_red    = [0.85, 0.33, 0.10];
color_yellow = [0.93, 0.69, 0.13];
color_palette = [color_red; color_blue; color_yellow];

tl = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% -------------------------------------------------------------------------
% Subplot (a): Full Augmented Spectrum & Cluster Identification
% -------------------------------------------------------------------------
ax1 = nexttile(1);
hold(ax1, "on"); grid(ax1, "on"); box(ax1, "on");

yline(ax1, 0, 'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off');

% Draw translucent patches for each identified spurious cluster
offset = 1.5; 
for j = 1:Q
    val_max = best_spurious_mat(1, j);
    val_min = best_spurious_mat(n_sys, j);
    
    idx_start = find(LEs == val_max, 1);
    idx_end   = find(LEs == val_min, 1, 'last');
    
    patch(ax1, [idx_start-0.4, idx_end+0.4, idx_end+0.4, idx_start-0.4], ...
              [val_min-offset, val_min-offset, val_max+offset, val_max+offset], ...
              [0.2 0.6 1], 'FaceAlpha', 0.2, 'EdgeColor', [0.2 0.6 1], ...
              'LineWidth', 1, 'HandleVisibility', 'off');
end

% Plot original computed nodes
h_raw = plot(ax1, (1:N)', LEs, 'o', 'MarkerSize', 5, 'MarkerFaceColor', color_blue, ...
    'MarkerEdgeColor', 'k', 'DisplayName', 'Computed full spectrum');

% Highlight the extracted isolated physical LEs
h_phys = plot(ax1, best_phys_idx(:), physical_LEs, 'p', 'MarkerSize', 10, ...
    'MarkerFaceColor', color_yellow, 'MarkerEdgeColor', 'k', 'LineWidth', 1.2, ...
    'DisplayName', 'Physical LEs');

% Dummy handle for the legend of the patches
h_box = patch(ax1, 'XData', NaN, 'YData', NaN, 'FaceColor', [0.2 0.6 1], ...
    'FaceAlpha', 0.2, 'EdgeColor', [0.2 0.6 1], 'DisplayName', 'Spurious LE clusters');

xlabel(ax1, 'Index', 'Interpreter', 'latex'); 
ylabel(ax1, '$\lambda$', 'Interpreter', 'latex');
legend([h_raw, h_box, h_phys], 'Location', 'southwest', 'Interpreter', 'latex', 'FontSize', 11);
xlim(ax1, [0, N+1]); ylim(ax1, [-80, 20]);

text(ax1, 0.03, 0.94, '(a)', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 16);

% -------------------------------------------------------------------------
% Subplot (b): Convergence History of Physical LEs
% -------------------------------------------------------------------------
ax2 = nexttile(2);
hold(ax2, "on"); grid(ax2, "on"); box(ax2, "on");

for i = 1:length(physical_LEs)
    plot(ax2, Time, LEs_history(:, best_phys_idx(i)), 'LineWidth', 1.5, 'Color', color_palette(i, :));
end

yline(ax2, 0, 'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
xlabel(ax2, 'Time (s)', 'Interpreter', 'latex'); 
ylabel(ax2, '$\lambda$', 'Interpreter', 'latex');
legend('$\lambda_1$', '$\lambda_2$', '$\lambda_3$', 'Interpreter', 'latex', 'Location', 'east');
legend(ax2, 'boxoff');
ylim(ax2, [-20, 5]);

text(ax2, 0.03, 0.94, '(b)', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 16);

% -------------------------------------------------------------------------
% Subplot (c): Experimental LLE Validation
% -------------------------------------------------------------------------
ax3 = nexttile(3);
hold(ax3, "on"); grid(ax3, "on"); box(ax3, "on");

% Plot the actual logarithmic distance evolution
plot(ax3, t_exp, log_dist, 'Color', color_red, 'LineWidth', 1.2); 

% Plot the linear regression fit
t_fit_line = t_exp(idx_fit);
val_fit_line = polyval(P_fit, t_fit_line);
plot(ax3, t_fit_line, val_fit_line, '--', 'Color', color_blue, 'LineWidth', 2);

% Annotate the estimated slope
text(ax3, 0.15, 0.35, sprintf('Slope $\\approx$ %.3f', LLE_value), ...
    'Units', 'normalized', 'Interpreter', 'latex', 'BackgroundColor', 'w', 'EdgeColor', 'k');

xlabel(ax3, 'Time (s)', 'Interpreter', 'latex');
ylabel(ax3, '$\ln || \mathbf{x}(t) - \bar{\mathbf{x}}(t) ||$', 'Interpreter', 'latex');

text(ax3, 0.03, 0.94, '(c)', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 16);

% -------------------------------------------------------------------------
% Subplot (d): 0-1 Test for Chaos (p-q Translation)
% -------------------------------------------------------------------------
ax4 = nexttile(4);
hold(ax4, "on"); grid(ax4, "on"); box(ax4, "on");

% Plot Brownian-like translation trajectory
plot(ax4, p_c(10:end), q_c(10:end), 'Color', color_blue, 'LineWidth', 0.8);

xlabel(ax4, '$p_c$', 'Interpreter', 'latex');
ylabel(ax4, '$q_c$', 'Interpreter', 'latex');

% Annotate the K value
text(ax4, 0.75, 0.10, sprintf('$K \\approx %.3f$', K_val), ...
    'Units', 'normalized', 'Interpreter', 'latex', 'BackgroundColor', 'w', 'EdgeColor', 'k');

axis(ax4, 'equal');  % Crucial for correct visual assessment of diffusion
text(ax4, 0.03, 0.94, '(d)', 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 16);

% =========================================================================
% Export Graphic
% =========================================================================
exportgraphics(gcf, './Lorenz_ISR_Results_alpha_0.97.png', 'Resolution', 600);