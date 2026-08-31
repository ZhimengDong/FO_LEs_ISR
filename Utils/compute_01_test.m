function [K_median, p_selected, q_selected] = compute_01_test(phi)
%COMPUTE_01_TEST Apply the correlation form of the 0-1 test for chaos.

    phi = phi(:)';

    downsample_factor = 10;
    phi = phi(1:downsample_factor:end);
    phi = phi - mean(phi);

    sample_count = numel(phi);
    n_max = floor(sample_count / 10);
    if n_max < 2
        error('The supplied time series is too short for the 0-1 test.');
    end

    c_values = linspace(pi / 5, 4 * pi / 5, 100);
    K_values = zeros(size(c_values));
    n_values = 1:n_max;
    p_selected = [];
    q_selected = [];

    for i = 1:numel(c_values)
        c = c_values(i);
        indices = 1:sample_count;
        p = cumsum(phi .* cos(indices * c));
        q = cumsum(phi .* sin(indices * c));

        if i == round(numel(c_values) / 2)
            p_selected = p;
            q_selected = q;
        end

        mean_square_displacement = zeros(1, n_max);
        for n = 1:n_max
            mean_square_displacement(n) = mean( ...
                (p(1+n:end) - p(1:end-n)).^2 + ...
                (q(1+n:end) - q(1:end-n)).^2);
        end

        correlation = corrcoef(n_values, mean_square_displacement);
        K_values(i) = correlation(1, 2);
    end

    K_median = median(K_values);
end

