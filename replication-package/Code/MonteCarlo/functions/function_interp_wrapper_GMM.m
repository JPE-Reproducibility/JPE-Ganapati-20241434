function p = function_interp_wrapper_GMM(o, d, lnn, do_GFT)
  
% theta_true
    % p.theta_n = d.theta_true;
    p.theta_n = function_value_interp(d.R_nE_ij, o.theta, lnn);
    p.rho_n = function_value_interp(d.R_nE_ij, o.rho_el, lnn);
    p.epsilon_n = function_value_interp(d.R_nE_ij, o.epsilon_el, lnn);

    X_ij   = zeros(max(d.I), max(d.J)); idx  = sub2ind(size(X_ij), d.I, d.J); X_ij(idx) = exp(d.X_ij);
    n_ij   = zeros(max(d.I), max(d.J)); idx  = sub2ind(size(n_ij), d.I, d.J); n_ij(idx) = exp(d.R_nE_ij);
    
%     % define elasticity functions (estimated)
%     lnn         = @(n) evknots_make(d.k,log(n)) ;
%     rho         = @(n) exp((lnn(n)*o.est_rho')');                             
%     epsilon     = @(n) exp((lnn(n)*o.est_epsilon')');

    if do_GFT == 1
        % GFT (estimated)
        tic
        p.w_hat = GFT_fct_general_fast(o.epsilon_est_f , o.rho_est_f, d.sigma, X_ij, n_ij);
        toc
    
        % call elasticity functions (true)
        tic
        p.w_hat_true = GFT_fct_general_fast(d.epsilon_true.Value, d.rho_true.Value, d.sigma, X_ij, n_ij);
        toc
        % p.w_hat_true = d.w_hat_true;
    
        % GFT (true)
        p.w_hat_mse = mean((p.w_hat - p.w_hat_true).^2);
    end

    p.theta_true = function_value_interp(exp(d.R_nE_ij), d.theta_true, exp(lnn));
    % p.rho_true = function_value_interp(exp(d.R_nE_ij), d.rho_true, exp(lnn));
    % p.epsilon_true = function_value_interp(exp(d.R_nE_ij), d.epsilon_true, exp(lnn));
    p.rho_true_el = function_value_interp(exp(d.R_nE_ij), d.rho_true_el, exp(lnn));
    p.epsilon_true_el = function_value_interp(exp(d.R_nE_ij), d.epsilon_true_el, exp(lnn));
end