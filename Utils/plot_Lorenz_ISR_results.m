function [physical_LEs, K_value] = plot_Lorenz_ISR_results( ...
    result_path, output_path, delta)
%PLOT_LORENZ_ISR_RESULTS Recreate the four-panel fractional Lorenz figure.

    if nargin < 1 || isempty(result_path)
        result_path = fullfile( ...
            'Data', 'Lorenz_ISR_Results_alpha_0.995.mat');
    end
    if nargin < 2 || isempty(output_path)
        output_path = 'Lorenz_ISR_Results.png';
    end
    if nargin < 3 || isempty(delta)
        delta = 1e-1;
    end

    loaded = load(result_path);
    if ~isfield(loaded, 'Lorenz_ISR_Results')
        error('Expected variable Lorenz_ISR_Results in %s.', result_path);
    end
    results = loaded.Lorenz_ISR_Results;

    X_series = results.X_seq;
    time_axis = results.Time(:);
    LE_history = results.LEs_history;
    t_experimental = results.t_exp(:);
    log_distance = results.log_dist(:);
    LLE_value = results.LLE_value;
    fit_indices = results.idx_fit;
    fit_coeff = results.P_fit;

    [sorted_augmented, sort_order] = sort( ...
        LE_history(end, :), 2, 'descend');
    sorted_augmented = sorted_augmented(:);
    sort_order = sort_order(:);
    [physical_LEs, physical_indices, spurious_clusters, diagnostics] = ...
        extract_physical_LEs(sorted_augmented, 3, delta);
    physical_history_columns = sort_order(physical_indices);

    [K_value, p_translation, q_translation] = compute_01_test(X_series);

    fprintf('Data file: %s\n', result_path);
    fprintf('Degeneracy threshold delta: %.3g\n', delta);
    fprintf('Physical LEs: [%.6f, %.6f, %.6f]\n', physical_LEs);
    fprintf('Two-trajectory LLE: %.6f\n', LLE_value);
    fprintf('0-1 test K-median: %.6f\n', K_value);

    figure_handle = figure( ...
        'Position', [50, 50, 900, 750], ...
        'Color', 'w');
    set(0, 'DefaultAxesFontName', 'Times New Roman');
    set(0, 'DefaultAxesFontSize', 12);
    set(0, 'DefaultTextFontName', 'Times New Roman');
    set(0, 'DefaultTextFontSize', 12);

    blue = [0.00, 0.45, 0.74];
    red = [0.85, 0.33, 0.10];
    yellow = [0.93, 0.69, 0.13];
    colors = [red; blue; yellow];

    tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    %% (a) Augmented spectrum and extracted physical exponents
    ax1 = nexttile(1);
    hold(ax1, 'on');
    grid(ax1, 'on');
    box(ax1, 'on');
    yline(ax1, 0, 'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off');

    patch_offset = 1.5;
    for j = 1:size(spurious_clusters, 2)
        cluster_max = max(spurious_clusters(:, j));
        cluster_min = min(spurious_clusters(:, j));
        [~, start_index] = min(abs(sorted_augmented - cluster_max));
        [~, end_index] = min(abs(sorted_augmented - cluster_min));
        patch(ax1, ...
            [start_index-0.4, end_index+0.4, end_index+0.4, start_index-0.4], ...
            [cluster_min-patch_offset, cluster_min-patch_offset, ...
             cluster_max+patch_offset, cluster_max+patch_offset], ...
            [0.2, 0.6, 1.0], ...
            'FaceAlpha', 0.2, ...
            'EdgeColor', [0.2, 0.6, 1.0], ...
            'LineWidth', 1, ...
            'HandleVisibility', 'off');
    end

    raw_handle = plot(ax1, 1:numel(sorted_augmented), sorted_augmented, ...
        'o', 'MarkerSize', 5, 'MarkerFaceColor', blue, ...
        'MarkerEdgeColor', 'k', 'DisplayName', 'Computed full spectrum');
    physical_handle = plot(ax1, physical_indices, physical_LEs, ...
        'p', 'MarkerSize', 10, 'MarkerFaceColor', yellow, ...
        'MarkerEdgeColor', 'k', 'LineWidth', 1.2, ...
        'DisplayName', 'Physical LEs');
    cluster_handle = patch(ax1, ...
        'XData', NaN, 'YData', NaN, ...
        'FaceColor', [0.2, 0.6, 1.0], ...
        'FaceAlpha', 0.2, ...
        'EdgeColor', [0.2, 0.6, 1.0], ...
        'DisplayName', 'Spurious LE clusters');

    xlabel(ax1, 'Index');
    ylabel(ax1, '$\lambda$', 'Interpreter', 'latex');
    legend(ax1, [raw_handle, cluster_handle, physical_handle], ...
        'Location', 'southwest', 'Interpreter', 'latex', 'FontSize', 10);
    xlim(ax1, [0, numel(sorted_augmented) + 1]);
    ylim(ax1, [-80, 20]);
    text(ax1, 0.03, 0.94, '(a)', 'Units', 'normalized', ...
        'FontWeight', 'bold', 'FontSize', 16);

    %% (b) Convergence of the physical exponents
    ax2 = nexttile(2);
    hold(ax2, 'on');
    grid(ax2, 'on');
    box(ax2, 'on');
    for i = 1:numel(physical_LEs)
        plot(ax2, time_axis, LE_history(:, physical_history_columns(i)), ...
            'LineWidth', 1.5, 'Color', colors(i, :));
    end
    yline(ax2, 0, 'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
    xlabel(ax2, 'Time');
    ylabel(ax2, '$\lambda$', 'Interpreter', 'latex');
    legend(ax2, '$\lambda_1$', '$\lambda_2$', '$\lambda_3$', ...
        'Interpreter', 'latex', 'Location', 'east');
    ylim(ax2, [-20, 5]);
    text(ax2, 0.03, 0.94, '(b)', 'Units', 'normalized', ...
        'FontWeight', 'bold', 'FontSize', 16);

    %% (c) Two-trajectory LLE check
    ax3 = nexttile(3);
    hold(ax3, 'on');
    grid(ax3, 'on');
    box(ax3, 'on');
    plot(ax3, t_experimental, log_distance, ...
        'Color', red, 'LineWidth', 1.2);
    fit_time = t_experimental(fit_indices);
    plot(ax3, fit_time, polyval(fit_coeff, fit_time), '--', ...
        'Color', blue, 'LineWidth', 2);
    text(ax3, 0.15, 0.35, sprintf('Slope $\\approx$ %.3f', LLE_value), ...
        'Units', 'normalized', 'Interpreter', 'latex', ...
        'BackgroundColor', 'w', 'EdgeColor', 'k');
    xlabel(ax3, 'Time');
    ylabel(ax3, '$\ln\|\mathbf{x}(t)-\bar{\mathbf{x}}(t)\|$', ...
        'Interpreter', 'latex');
    text(ax3, 0.03, 0.94, '(c)', 'Units', 'normalized', ...
        'FontWeight', 'bold', 'FontSize', 16);

    %% (d) 0-1 test
    ax4 = nexttile(4);
    hold(ax4, 'on');
    grid(ax4, 'on');
    box(ax4, 'on');
    plot(ax4, p_translation(10:end), q_translation(10:end), ...
        'Color', blue, 'LineWidth', 0.8);
    xlabel(ax4, '$p_c$', 'Interpreter', 'latex');
    ylabel(ax4, '$q_c$', 'Interpreter', 'latex');
    text(ax4, 0.72, 0.10, sprintf('$K \\approx %.3f$', K_value), ...
        'Units', 'normalized', 'Interpreter', 'latex', ...
        'BackgroundColor', 'w', 'EdgeColor', 'k');
    axis(ax4, 'equal');
    text(ax4, 0.03, 0.94, '(d)', 'Units', 'normalized', ...
        'FontWeight', 'bold', 'FontSize', 16);

    exportgraphics(figure_handle, output_path, 'Resolution', 600);
    fprintf('Figure written to: %s\n', output_path);
    fprintf('Extraction cost: %.9g\n', diagnostics.best_cost);
end
