function y = logpareto_inv_shifted(p, k, alpha, mu)
% LOGPARETO_INV_SHIFTED computes the inverse CDF of the shifted log-Pareto.
%
%   y = LOGPARETO_INV_SHIFTED(p, k, alpha, mu) returns the inverse CDF 
%   of the log-Pareto distribution with scale k, shape alpha, and 
%   location mu, evaluated at the probabilities in p.

    if any(p < 0 | p > 1)
        error('Probabilities must be between 0 and 1.');
    end
    
%     if k <= 0 || alpha <= 0
%         error('Scale (k) and shape (alpha) parameters must be positive.');
%     end

    % Calculate the standard log-Pareto value
    x = exp(k ./ (1 - p).^(1/alpha));
    
    % Add the location parameter to shift the distribution
    y = x + mu;
end