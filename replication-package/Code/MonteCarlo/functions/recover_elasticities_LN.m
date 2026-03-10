function [epsilon, rho, epsilon_el, rho_el] = recover_elasticities_LN(o,d)
        
        n_sigma = o.est_sigma;
        
        %epsilon_n = @(n,nu) exp(nu .* norminv(1 - n, 0, 1));
        eps = 1e-3;
        %% TODO precompute integral
        %rho_grid = arrayfun(@(lb, ub) integral(epsilon_n, lb, ub), zeros(size(lnn, 1), 1), exp(lnn)); 
        %rho_n = @(n) interp1(exp(lnn), rho_grid, n, 'linear', 'extrap'); 
        z_f = @(n) exp(sqrt(2)*n_sigma*erfinv(2*(1-n)-1));
        e_f = @(n) (n);
        e = @(n) e_f(z_f(n));
        
        % Use integral with relaxed tolerances and better error handling
        r_f = @(n) arrayfun(@(n_val) compute_r_integral(n_val, z_f, e_f), n);
        
        function result = compute_r_integral(n_val, z_f, e_f)
            if n_val <= 0 || isnan(n_val)
                result = NaN;
                return;
            end
            
            % Integration limits - stay away from n=1 boundary
            n_lower = max(1e-10, n_val * 1e-8);
            n_upper = min(n_val, 1 - 1e-8);
            
            if n_upper <= n_lower
                result = NaN;
                return;
            end
            
            try
                % Relaxed tolerances for difficult integrals
                integral_val = integral(@(n2) e_f(z_f(n2)), n_lower, n_upper, ...
                    'RelTol', 1e-4, ...      % Relative tolerance (default 1e-6)
                    'AbsTol', 1e-8, ...      % Absolute tolerance (default 1e-10)
                    'ArrayValued', true);
                
                result = integral_val / n_val;
            catch ME
                % If integration fails, try with even more relaxed tolerances
                try
                    integral_val = integral(@(n2) e_f(z_f(n2)), n_lower, n_upper, ...
                        'RelTol', 1e-3, ...
                        'AbsTol', 1e-6, ...
                        'ArrayValued', true);
                    result = integral_val / n_val;
                catch
                    result = NaN;
                end
            end
        end
        
        r_ex = @(n) r_f(n);
        le = @(lnn) log(e(exp(lnn)));
        lr = @(lnn) log(r_ex(exp(lnn)));
        nvals = exp(d.R_nE_ij);
        
        epsilon_el = (le(log(nvals)+eps) - le(log(nvals)))./(eps);
        rho_el = (lr(log(nvals)+eps) - lr(log(nvals)))./(eps);
        epsilon = exp(le(log(nvals)));
        rho = exp(lr(log(nvals)));
end