function [omega_k, eta_k] = generate_oustaloup_params( ...
    alpha, recursion_order, w_low, w_high)
%GENERATE_OUSTALOUP_PARAMS Generate ISR nodes and spectral weights.

    if recursion_order < 1 || recursion_order ~= floor(recursion_order)
        error('recursion_order must be a positive integer.');
    end
    if w_low <= 0 || w_high <= w_low
        error('Require 0 < w_low < w_high.');
    end

    mu = (w_high / w_low)^((1 - alpha) / recursion_order);
    nu = (w_high / w_low)^(alpha / recursion_order);

    omega_zero = zeros(recursion_order, 1);
    omega_pole = zeros(recursion_order, 1);

    omega_zero(1) = w_low * sqrt(nu);
    omega_pole(1) = mu * omega_zero(1);

    for k = 1:recursion_order-1
        omega_zero(k+1) = nu * omega_pole(k);
        omega_pole(k+1) = mu * omega_zero(k+1);
    end

    gain = 10^((1 - alpha) * log10(w_low));
    residues = zeros(recursion_order, 1);

    for i = 1:recursion_order
        numerator = gain * ...
            (omega_pole(i) - omega_zero(i)) / omega_zero(i);
        product_term = 1;
        for j = 1:recursion_order
            if j ~= i
                product_term = product_term * ...
                    (omega_pole(i) - omega_zero(j)) / ...
                    (omega_pole(i) - omega_pole(j));
            end
        end
        residues(i) = numerator * product_term;
    end

    omega_k = [0; omega_pole(:)];
    eta_k = [gain; residues(:)];
end

