function [omega_k, c_k] = generate_oustaloup_params(alpha, N, w_low, w_high)
    % GENERATE_OUSTALOUP_PARAMS Generates nodes and weights for the ISR framework
    % using the Oustaloup recursive approximation method (Trigeassou variant).
    % (Method A in the associated manuscript).
    %
    % Inputs:
    %   alpha  : Fractional order
    %   N  : Number of recursive steps (Total nodes Q = N + 1)
    %   w_low  : Lower bound of the approximation frequency band
    %   w_high  : Upper bound of the approximation frequency band
    %
    % Outputs:
    %   omega_k  : Frequency nodes (poles of the integer-order proxy system)
    %   c_k  : Spectral weights (residues for the pseudo-state projection)
    
    mu = (w_high/w_low)^((1-alpha)/N);
    nu = (w_high/w_low)^(alpha/N);
    
    omega_prime = zeros(N, 1);
    omega_poles = zeros(N, 1);
    
    omega_prime(1) = w_low * sqrt(nu);
    omega_poles(1) = mu * omega_prime(1);
    
    for k = 1:N-1
        omega_prime(k+1) = nu * omega_poles(k);
        omega_poles(k+1) = mu * omega_prime(k+1);
    end
    
    % Global gain factor
    Gn = 10^((1-alpha)*log10(w_low)); 
    
    c_k = zeros(N, 1);
    % Calculate the residue for each pole using partial fraction expansion
    for i = 1:N
        num = Gn * (omega_poles(i) - omega_prime(i)) / omega_prime(i);
        prod_term = 1;
        for j = 1:N
            if j ~= i
                term = (1 - omega_poles(i)/omega_prime(j)) / (1 - omega_poles(i)/omega_poles(j));
                prod_term = prod_term * term;
            end
        end
        c_k(i) = num * prod_term;
    end
    
    % Append the integrating mode (omega = 0)
    omega_k = [0; omega_poles];
    c_k = [Gn; c_k];
    
    % Ensure outputs are strictly column vectors
    omega_k = omega_k(:);
    c_k = c_k(:);
end