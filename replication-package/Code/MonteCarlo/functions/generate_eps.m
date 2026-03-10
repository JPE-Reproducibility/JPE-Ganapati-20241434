function e = generate_eps(u, e_bar, alpha_e, gamma_e)
    % Solves G^e(e) = u for e
    % G^e(e) = 1 - (e/e_bar)^(-alpha_e) * (ln(e)/ln(e_bar))^(-gamma_e)
    % Define the equation to solve: G^e(e) - u = 0
    cdf_eq = @(e) 1 - (e/e_bar)^(-alpha_e) * (log(e)/log(e_bar))^(-gamma_e) - u;
    
    % Use fzero to find the root
    try
        e = fzero(cdf_eq, [e_bar, 5e7]);
    catch
        try
            e = fzero(cdf_eq, [e_bar, 1e3]);
        catch
            
            % warning('Could not solve for epsilon at u = %.6f', u);
            e = e_bar;
            %e = 1e3;

        end
    end
end