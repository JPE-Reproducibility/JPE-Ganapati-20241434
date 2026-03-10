function p = function_interp_wrapper(o, d, lnn, do_GFT)
  

    p.theta_n = function_value_interp(d.R_nE_ij, o.theta, lnn);
    p.rho_n = function_value_interp(d.R_nE_ij, o.rho_el, lnn);
    p.epsilon_n = function_value_interp(d.R_nE_ij, o.epsilon_el, lnn);
    %o_true = o;
    X_ij      = zeros(max(d.I), max(d.J));  
    n_ij      = zeros(max(d.I), max(d.J)); 
    idx       = sub2ind(size(X_ij), d.I, d.J);     X_ij(idx) = exp(d.X_ij);
    idx       = sub2ind(size(n_ij), d.I, d.J);     n_ij(idx) = exp(d.R_nE_ij);

    %% Test alt GFT if QQ Estimator
    n_mu        = o.ests2(1);
    n_sigma     = o.ests2(2); 

    ln_grid = linspace(-20, -.0001, 2000)'; 
    % extensive margin elasticity function
    z_f         = @(n) exp(sqrt(2)*n_sigma*erfinv(2*(1-n)-1));
    epsilon     = @(n) (z_f(exp(log(n))));
    r_f         = @(n) 1./n .* arrayfun(@(l,u) integral( @(n2) (z_f(n2)) ,l,u),zeros(size(n)),n);
    r_f_grid    = arrayfun(@(ln) r_f(exp(ln)), ln_grid);
    rho         = @(n) interp1(ln_grid, r_f_grid, log(n), 'spline', 'extrap'); % why are we interpolating. i need a good reason.


    if do_GFT == 1
        %% GFT (estimated)
        tic
        p.w_hat = GFT_fct_general_fast(epsilon, rho, d.sigma, X_ij, n_ij);
        toc
    
        %% GFT (truth)
        tic
        p.w_hat_true = GFT_fct_general_fast(d.epsilon_true.Value, d.rho_true.Value, d.sigma, X_ij, n_ij);
        % p.w_hat_true = GFT_fct_general_fast(epsilonT, rhoT, d.sigma, X_ij, n_ij);
        toc
       
        % scatter(p.w_hat,p.w_hat_true)
        % GFT (true)
        p.w_hat_mse = mean((p.w_hat - p.w_hat_true).^2);
        disp("GFT done")

    end

    p.theta_true = function_value_interp(exp(d.R_nE_ij), d.theta_true, exp(lnn));
    % p.rho_true = function_value_interp(exp(d.R_nE_ij), d.rho_true, exp(lnn));
    % p.epsilon_true = function_value_interp(exp(d.R_nE_ij), d.epsilon_true, exp(lnn));
    p.rho_true_el = function_value_interp(exp(d.R_nE_ij), d.rho_true_el, exp(lnn));
    p.epsilon_true_el = function_value_interp(exp(d.R_nE_ij), d.epsilon_true_el, exp(lnn));
end