function solve_FP(project_path, output_path, data_params, epsilon, eps_true, rho, r_true, kappa, H_e, Q_lnx, p_grid, sims, cores, upd_step)
    % Extract parameters from struct
    fct_form = data_params.fct_form;
    sigma = data_params.sigma;
    countries = data_params.c;
    mu_f = data_params.mu_f;
    sigma_z = data_params.sigma_z;
    sigma_r = data_params.sigma_r;
    sigma_f = data_params.sigma_f;
    kappa_r = data_params.kappa_r;
    kappa_f = data_params.kappa_f;
    kappa_r_home = data_params.kappa_r_home;
    kappa_f_home = data_params.kappa_f_home;
    mu_r = data_params.mu_r;
    mu_z = data_params.mu_z;
    addpath(project_path + '/functions', '-frozen')
    data_path = output_path + "/data";
    
    % Countries
    c = countries;
    % # of simulations
    J = sims;
    
    Lbar = ones(c, 1);       % Labor endowment
    Fbar = ones(c, 1);       % 
    
    tol = 1e-7;
    error_p = 100;
    %error_w = 100;
    gamma = upd_step; %contraction mapping
    
    seeds        = (123+(1:J))';
    conv_flag    = false(J,1);
    nan_flag     = false(J,1);
    last_err_p   = nan(J,1);
    last_err_w   = nan(J,1);
    err_flag     = false(J,1);
    err_id       = strings(J,1);
    err_msg      = strings(J,1);
    
    %% Equilibrium objects
    n_ij        = @(fbar, rbar, w_i, w_j, P_j) max(min(1 - H_e(log((sigma .* fbar .* w_i.^sigma)./(rbar .* P_j .* w_j))), 0.999), 1e-18);
    xbar_ij     = @(rbar, rho, w_i, w_j, P_j) (rbar .* rho .* w_i.^(1 - sigma) .* P_j .* w_j);
    N_i         = @(fbarmat, kappa) Lbar ./ (sigma * (Fbar + sum(fbarmat .* kappa, 2)));

    %% Initialization
    R_xbar_ij = zeros(c^2, J); % log average revenue
    R_n_ij = zeros(c^2, J); % log probability of entry
    Z1 = zeros(c^2, J); % bilateral shifter
    Z2 = zeros(c^2, J); % bilateral shifter 2 (=1)
    X_ij = zeros(c^2, J); % trade flows
    w_i_sim = zeros(c, J); % wages
    P_i_sim = zeros(c, J); % prices
    N_i_sim = zeros(c, J); % firm entry/mass of new firms
    GFT_T = ones(c, J); % GFT (unused for now)
    rho_true = zeros(c^2, J); % equilibrium elasticity function values
    epsilon_true = zeros(c^2, J); % equilibrium elasticity function values
    theta_true = zeros(c^2, J); % trade elasticity
    rho_true_el = zeros(c^2, J); % intensive margin elasticity
    epsilon_true_el = zeros(c^2, J); % extensive margin elasticity
    
    pp = numel(p_grid); % number of quantiles of sales distribution
    Q_lnx_all = nan(c^2, pp, J); % sales
    
    [i_idx, j_idx] = ndgrid(1:c, 1:c); 
    i_idx = repmat(i_idx(:), 1, J);  % c^2 x 1 %origin id
    j_idx = repmat(j_idx(:), 1, J);  % c^2 x 1 %destination id
    home_idx = (1:c+1:c^2)';

    delete(gcp('nocreate'))
    parpool(min(cores, sims));
    tic;
    parfor b = 1:J
        %% Baseline equilibrium
      try
        % reproducibility under parallelization
        rng(seeds(b));
        
        z = mu_z + randn(c, c).*sigma_z; %bilateral revenue shifter
       
        eps_r = randn(c, c).*sigma_r + mu_r; %rbar error term              
        eps_f = unifrnd(mu_f - sigma_f, mu_f + sigma_f, c); %fbar error term
        % revenue potentials
        rbar = exp(kappa_r*z + kappa_r_home*eye(c) + eps_r);
        rbar = rbar(:);           % flatten to c^2 x 1
        % fixed costs
        fbar = max(1e-7, eps_f.*exp(kappa_f*z - kappa_f_home*eye(c)));
        fbar = fbar(:);           % c^2 x 1
    
        z = z(:);
        
        w = ones(c, 1);            % initial wages
        P = ones(c, 1);            % initial prices 
        i_outer = 0;
        i_inner = 0;
        error_w = 100;
       
        while error_w > tol
            error_p = 100;
            w_i_vec = kron(ones(c,1), w);                           % c^2 x 1
            w_j_vec   = kron(w, ones(c,1));                         % c^2 x 1
            while error_p > tol
                % prices at dest
                
                P_j_vec   = kron(P.^(sigma - 1), ones(c,1));            % c^2 x 1
                
                % probabilities of entry n_{ij}
                nvec = n_ij(fbar, rbar, w_i_vec, w_j_vec, P_j_vec);
                nmat = reshape(nvec, c, c)
                % intesive margin elasticity function
                rho_vec = rho(nvec);
                % average firm sales xbar_{ij}
                xbar_vec = xbar_ij(rbar, rho_vec, w_i_vec, w_j_vec, P_j_vec);
               
            
                % N_i
                fbarmat = reshape(fbar, c, c);                          % c x c
                kappa_vec = kappa(nvec);
                kappa_mat = reshape(kappa_vec, c, c);
                N = N_i(fbarmat, kappa_mat);
                %fprintf('N_i: %.3e\n', any(isnan(N)));
            
                % total good demand \sum_i N_i n_{ij} xbar_{ij}
                N_rep = kron(ones(c,1), N);                             % c^2 x 1
                Xvec = N_rep .* nvec .* xbar_vec;                      % c^2 x 1
                Xmat = reshape(Xvec, c, c);                           % c x c
            
                tgd = sum(Xmat, 1)';                                   % c x 1
            
                % excess good demand
                egd = w - tgd;                                          % c x 1
            
                % update prices
                P_new = P + gamma * egd;
            
                
                error_p = max(abs(P_new - P));
                P = P_new;
            
                %fprintf('Error: %.3e\n', error_p);
                i_inner = i_inner + 1;
            end
        
        % total labor demand
        tld = sum(Xmat, 2);
        
        % excess labor demand
        eld = tld - w; 
        %wlog normalize relative to the first country
        eld = eld - eld(1);
        
        % update wages
        w_new = w + gamma*eld;
        
        error_w = max(abs(w_new - w));
        w = w_new;
        w(1) = 1;
            
        fprintf('Error: %.3e\n', error_w);
        
        i_outer = i_outer + 1;
        
        end

        fprintf('Converged');

        %% non-converged seed values
        conv_flag(b)  = (error_w <= tol);   % inner loop already broke when error_p <= tol
        last_err_p(b) = error_p;
        last_err_w(b) = error_w;
        % identify any non-finite along the key vectors/matrices
        bad_any = ~all(isfinite([nvec; xbar_vec; Xvec; N; P; w]));
        nan_flag(b) = bad_any;
        
        if bad_any
           badX = ~isfinite(Xvec); badN = ~isfinite(N); badP = ~isfinite(P); badw = ~isfinite(w);
           fprintf('  offenders: X=%d N=%d P=%d w=%d\n', any(badX), any(badN), any(badP), any(badw));
        end
        
        %% Collect objects
        R_n_ij(:, b) = log(nvec);
        R_xbar_ij(:, b) = log(xbar_vec);
        X_ij(:, b) = log(Xvec);
        Z1(:, b) = z;
        Z2(:, b) = z;

        rho_true_el(:, b) = r_true(nvec);
        epsilon_true_el(:, b) = eps_true(nvec);
        theta_true(:, b) = (sigma - 1).*(1 - ((1 + rho_true_el(:, b))./epsilon_true_el(:, b)));
        
        Rbar_vec = xbar_vec ./ rho_vec;
        Q_lnx_all(:, :, b) = Q_lnx(Rbar_vec, nvec);
    
        P_i_sim(:, b) = P;
        w_i_sim(:, b) = w;
        N_i_sim(:, b) = N;
        
        
        catch ME
        % mark run as an error 
        err_flag(b) = true;
        err_id(b)   = string(ME.identifier);
        err_msg(b)  = string(ME.message);
        conv_flag(b)= false;
        nan_flag(b) = true;  
        fprintf('[run %d] ERROR: %s - %s\n', b, ME.identifier, ME.message);
        continue
        
      end
    end
    toc;
    
    keep = conv_flag & ~nan_flag & ~err_flag;
    
    % objects for M
    R_n_ij    = R_n_ij(:, keep);
    R_xbar_ij = R_xbar_ij(:, keep);
    Z1        = Z1(:, keep);
    Z2        = Z2(:, keep);
    X_ij      = X_ij(:, keep);
    Q_lnx_all = reshape(Q_lnx_all, c^2.*pp, J);
    Q_lnx_all = Q_lnx_all(:, keep);
    
    
    i_idx_keep = i_idx(:, 1:sum(keep));   % or just reuse i_idx(:,keep) 
    j_idx_keep = j_idx(:, 1:sum(keep));
    
    %GFT_T   = GFT_T(:, keep);
    P_i_sim = P_i_sim(:, keep);
    w_i_sim = w_i_sim(:, keep);
    N_i_sim = N_i_sim(:, keep);
    
    %if exist('rho_true','var'),     rho_true     = rho_true(:, keep);     end
    %if exist('epsilon_true','var'), epsilon_true = epsilon_true(:, keep); end
    rho_true = @(n) rho(n);
    epsilon_true = @(n) epsilon(log(n));
    
    if exist('theta_true','var'),   theta_true   = theta_true(:, keep);   end
    if exist('epsilon_true_el','var'),   epsilon_true_el   = epsilon_true_el(:, keep);   end
    if exist('rho_true_el','var'),   rho_true_el   = rho_true_el(:, keep);   end
    
    if any(keep)
        M = cat(3, i_idx_keep, j_idx_keep, R_n_ij, R_xbar_ij, Z1, Z2, X_ij);
        M = permute(M, [1, 3, 2]);   % (c^2, K, #kept)
    else
        M = zeros(c^2, 7, 0);        % empty placeholder if nothing kept
    end
    
    % ---- run summary - identify converged/non-converged LN
    run_summary = table( ...
    (1:J)', seeds, conv_flag, nan_flag, err_flag, last_err_p, last_err_w, err_id, err_msg, ...
    'VariableNames', {'iter','seed','converged','hasNaN','hadError','last_err_p','last_err_w','error_id','error_msg'});
    
    % store params
    q.sigma_true     = sigma;
    q.sigma_z        = sigma_z;
    q.sigma_r        = sigma_r;
    q.sigma_f        = sigma_f;
    q.kappa_r        = kappa_r;
    q.kappa_r_home   = kappa_r_home;
    q.kappa_f        = kappa_f;
    q.kappa_f_home   = kappa_f_home;
    q.mu_f           = mu_f;

    % ---- save runs ----
   save(fullfile(data_path, "Simulated_data_" + fct_form + "_alt.mat"), ...
         "M","N_i_sim","P_i_sim","w_i_sim", ...
         "rho_true","epsilon_true","rho_true_el","epsilon_true_el","theta_true","q", "Q_lnx_all", "pp");
    save(fullfile(data_path, "Run_sum" + fct_form + "_alt.mat"), "run_summary");

end