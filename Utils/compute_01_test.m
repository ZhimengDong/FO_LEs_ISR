function [K_median, p_qc, q_qc] = compute_01_test(phi)
    % COMPUTE_01_TEST Evaluates the 0-1 test for chaos on a 1D time series.
    % Returns the asymptotic correlation coefficient K near 1 (chaos) or 0 (regular).
    
    phi = phi(:)'; % Ensure row vector
    
    % Downsampling: Adjusts the sampling interval closer to ~1/10 of the characteristic period
    ds_factor = 10; 
    phi_ds = phi(1:ds_factor:end);
    
    % Mean Subtraction
    % Removing the mean forces the complex oscillatory integral term to vanish,
    phi_ds = phi_ds - mean(phi_ds);
    
    N = length(phi_ds);
    
    % Select translation constants c within [pi/5, 4pi/5] to avoid resonance
    N_c = 100;
    c_vals = linspace(pi/5, 4*pi/5, N_c); 
    K_vals = zeros(1, N_c);
    
    n_max = round(N / 10);
    n_arr = 1:n_max;
    
    p_qc = []; q_qc = [];
    
    for i = 1:N_c
        c = c_vals(i);
        
        % 1. Calculate translation variables p and q
        p = cumsum(phi_ds .* cos((1:N)*c));
        q = cumsum(phi_ds .* sin((1:N)*c));
        
        % Extract a specific (p,q) trajectory for phase portrait visualization (median c)
        if i == round(N_c / 2)
            p_qc = p; q_qc = q;
        end
        
        % 2. Calculate mean square displacement M_c
        M_c = zeros(1, n_max);
        for n = 1:n_max
            % Vectorized calculation of the mean Euclidean distance squared
            M_c(n) = mean((p(1+n:N) - p(1:N-n)).^2 + (q(1+n:N) - q(1:N-n)).^2);
        end
        
        % 3. Calculate K using the correlation coefficient method
        % (Since phi_ds is zero-mean, M_c is equivalent to the modified D_c)
        corr_mat = corrcoef(n_arr, M_c);
        K_vals(i) = corr_mat(1, 2);
    end
    
    % Return the median to eliminate outliers from rare resonances
    K_median = median(K_vals);
end