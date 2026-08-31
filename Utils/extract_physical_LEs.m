function [physical_LEs, physical_indices, spurious_clusters, diagnostics] = ...
    extract_physical_LEs(augmented_LEs, physical_dim, delta)
%EXTRACT_PHYSICAL_LES Isolate physical exponents by degeneracy rejection.
%   DELTA is the minimum normalization threshold in Eq. (31) of the paper.

    if nargin < 3 || isempty(delta)
        delta = 1e-1;
    end
    if delta <= 0
        error('delta must be positive.');
    end

    sorted_LEs = sort(augmented_LEs(:), 'descend');
    spectrum_size = numel(sorted_LEs);
    cluster_count = (spectrum_size - physical_dim) / physical_dim;

    if cluster_count ~= floor(cluster_count) || cluster_count < 1
        error('Spectrum size is incompatible with the physical dimension.');
    end

    candidates = nchoosek(1:spectrum_size, physical_dim);
    costs = inf(size(candidates, 1), 1);
    best_cost = inf;
    physical_indices = [];
    spurious_clusters = [];

    for i = 1:size(candidates, 1)
        candidate_indices = candidates(i, :);
        mask = true(spectrum_size, 1);
        mask(candidate_indices) = false;
        candidate_spurious = reshape( ...
            sorted_LEs(mask), physical_dim, cluster_count);

        cluster_width = max(candidate_spurious) - min(candidate_spurious);
        cluster_magnitude = max(abs(candidate_spurious));
        cost = sum(cluster_width ./ max(cluster_magnitude, delta));
        costs(i) = cost;

        if cost < best_cost
            best_cost = cost;
            physical_indices = candidate_indices;
            spurious_clusters = candidate_spurious;
        end
    end

    physical_LEs = sorted_LEs(physical_indices);
    diagnostics = struct();
    diagnostics.delta = delta;
    diagnostics.best_cost = best_cost;
    diagnostics.sorted_augmented_LEs = sorted_LEs;
    diagnostics.candidate_costs = costs;
end

