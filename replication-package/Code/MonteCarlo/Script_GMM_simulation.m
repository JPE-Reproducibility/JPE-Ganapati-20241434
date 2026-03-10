function Script_GMM_simulation(project_path, output_path, fct_form, estimators, cores, sims, do_GFT, varargin)
    % Parse optional arguments
    p = inputParser;
    addParameter(p, 'RunSingle', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x > 0));
    addParameter(p, 'RunRange', [], @(x) isempty(x) || (isnumeric(x) && numel(x) == 2));
    parse(p, varargin{:});
    
    run_single = p.Results.RunSingle;
    run_range = p.Results.RunRange;
    
    % Validate inputs
    if ~isempty(run_single) && ~isempty(run_range)
        error('Cannot specify both RunSingle and RunRange');
    end
    
    addpath(project_path + '/functions/', '-frozen')
    addpath(project_path + '/functions/GFT_v2', '-frozen')
    
    data_path = output_path + "/data";
    estimation_path = output_path + "/estimation";
    figure_path = output_path + "/plots";
    
    %% Load And Run
    
    Ncores = cores;
    Nsim = sims;
    
    % Determine which simulations to run
    if ~isempty(run_single)
        sim_indices = run_single;
        if run_single > Nsim
            error('RunSingle index %d exceeds total simulations %d', run_single, Nsim);
        end
        fprintf('Running SINGLE economy #%d on %d core(s)\n', run_single, Ncores);
    elseif ~isempty(run_range)
        sim_indices = run_range(1):run_range(2);
        if run_range(2) > Nsim
            error('RunRange upper bound %d exceeds total simulations %d', run_range(2), Nsim);
        end
        fprintf('Running economies %d to %d (%d total) on %d core(s)\n', ...
            run_range(1), run_range(2), numel(sim_indices), Ncores);
    else
        sim_indices = 1:Nsim;
        fprintf('Running on %d cores, %d total economies\n', Ncores, Nsim);
    end
    
    n_to_run = numel(sim_indices);
    max_index = max(sim_indices);  % Get the largest index we'll need
    
    % Load data
    simfile = fullfile(data_path, "Simulated_data_" + fct_form + "_alt.mat");
    load(simfile)   % expects M, epsilon_true, rho_true, theta_true
    
    %% Fixed params
    
    sigma     = 3.2;
    kappa_tau = 1;
    kappa_f   = 0;
    kappa_eps = 1/((sigma-1)*kappa_tau+kappa_f);

    filter     = find(ones(size(M(:,1)))==1);
    R_nE_ij    = squeeze((M(filter,3, :)));
    lnn        = linspace(prctile(R_nE_ij, 1, 'all'), prctile(R_nE_ij, 99, 'all'), 500)';

    for ee = 1:numel(estimators)
        basis = estimators{ee}.basis;
        knots = estimators{ee}.knots;

        fprintf('\n=== Running: fct_form=%s | basis=%s | knots=%d ===\n', fct_form, basis, knots);
        tic;
        
        % =========================
        % Parpool machinery – START
        
        % Broadcast big inputs once per worker
        M_const           = parallel.pool.Constant(@() M);
        eps_true_el_const = parallel.pool.Constant(@() epsilon_true_el);
        rho_true_el_const = parallel.pool.Constant(@() rho_true_el);
        theta_true_const  = parallel.pool.Constant(@() theta_true);
        % epsilon_true_const  = parallel.pool.Constant(@() epsilon_true);  %% SG added- i think needed for parllel?
        % rho_true_const  = parallel.pool.Constant(@() rho_true);  %% SG added- i think needed for parllel?
        eps_true_const    = parallel.pool.Constant(@() epsilon_true);
        rho_true_const    = parallel.pool.Constant(@() rho_true);

        if strcmp(fct_form,"Estimates")     
            baselineConst = parallel.pool.Constant(@() load(fullfile(data_path,'baseline_estimate.mat')));
        end

        % Decide whether to use parallel pool based on NUMBER OF CORES
        use_parallel = (Ncores > 1);
        
        if use_parallel
            % delete(gcp('nocreate'))
            % parpool(min(Ncores, n_to_run))
            currentPool = gcp('nocreate');
            if isempty(currentPool)
                parpool(min(Ncores, n_to_run))
            else
            end

            pctRunOnAll maxNumCompThreads(1);
            pctRunOnAll warning('off','MATLAB:integral:NonFiniteValue')
        else
            % For single core, just suppress warnings
            warning('off','MATLAB:integral:NonFiniteValue')
            fprintf('  Using sequential processing (1 core)\n');
        end
        % Parpool machinery – END
        % =========================

        % Preallocate cells based on the maximum index we'll access
        oC = cell(1, max_index);
        pC = cell(1, max_index);
        kC = cell(1, max_index);
        RnEijC = cell(1, max_index);
        converg = false(1, max_index);
        err = cell(max_index, 1);
        
        % Main loop - use parfor if Ncores > 1, regular for if Ncores == 1
        if use_parallel
            parfor i = sim_indices  % Iterate directly over sim_indices
                [oC{i}, pC{i}, kC{i}, RnEijC{i}, converg(i), err{i}] = ...
                    run_single_simulation(i, M_const, eps_true_const, rho_true_const, ...
                    eps_true_el_const, rho_true_el_const, theta_true_const, ...
                    knots, sigma, basis, lnn, kappa_tau, kappa_eps,do_GFT);
            end
        else
            for idx = 1:n_to_run
                i = sim_indices(idx);
                fprintf('  Processing economy %d/%d...\n', idx, n_to_run);
                [oC{i}, pC{i}, kC{i}, RnEijC{i}, converg(i), err{i}] = ...
                    run_single_simulation(i, M_const, eps_true_const, rho_true_const, ...
                    eps_true_el_const, rho_true_el_const, theta_true_const, ...
                    knots, sigma, basis, lnn,  kappa_tau, kappa_eps,do_GFT);
            end
        end
        disp("rows done")
        
        % Track converged/non-converged (only for the indices we ran)
        idx_converged = sim_indices(converg(sim_indices));
        idx_failed    = sim_indices(~converg(sim_indices));
        save(fullfile(estimation_path, "converged_idx.mat"), 'idx_converged','idx_failed','err');
        
        ok = ~cellfun(@isempty, pC);
        if ~any(ok)
            error('All iterations failed. First error:\n%s', err{find(~ok,1)});
        end
        ok_idx = find(ok);

        % Unpack p
        p  = struct();
        fn = fieldnames(pC{ok_idx(1)});
        for k = 1:numel(fn)
            pieces = cellfun(@(s) s.(fn{k}), pC(ok_idx), 'UniformOutput', false);
            if isnumeric(pieces{1})
                pieces = cellfun(@(x) reshape(x,[],1), pieces, 'UniformOutput', false);
            end
            p.(fn{k}) = [pieces{:}];
        end

        % Unpack o
        o  = struct();
        fn = fieldnames(oC{ok_idx(1)});
        for j = 1:numel(fn)
            o.(fn{j}) = cellfun(@(s) s.(fn{j}), oC(ok_idx), 'UniformOutput', false);
        end

        % Gen only R_nE_ij of successful runs
        Rn_mat = cat(2, RnEijC{ok_idx});
        
        if knots == 1
            k_mat = [];
        else
            firstNonEmpty = find(~cellfun(@isempty, kC(ok_idx)), 1, 'first');
            assert(~isempty(firstNonEmpty), 'Expected non-empty d.k when knots ~= 1.');
            
            Lk = numel(kC{ok_idx(firstNonEmpty)});
            assert(all(cellfun(@(v) ~isempty(v) && numel(v) == Lk, kC(ok_idx))), ...
                'd.k length varies or is empty across sims when knots ~= 1.');
        
            k_mat = cell2mat(cellfun(@(v) v(:), kC(ok_idx), 'UniformOutput', false));
        end

        % Log errors
        bad = find(~ok & ~cellfun(@isempty, err));
        if ~isempty(bad)
            logf = fullfile(estimation_path, sprintf('errors_%s_%s_knots%d.txt', fct_form, basis, knots));
            fid = fopen(logf,'w');
            for t = 1:numel(bad)
                printf(fid, 'i = %d\n%s\n\n', bad(t), err{bad(t)});
            end
            fclose(fid);
        end

        % Save per-run results
        save_name_core = "Simulated_" + fct_form + "_Estimates_" + basis + "_" + num2str(knots);
        o.rho_est_f = {};
        o.epsilon_est_f={};
        save(fullfile(estimation_path, save_name_core + "_alt.mat"), "p", "o", "Rn_mat", "k_mat","ok_idx","-v7.3");

        fprintf('Completed in %.2f seconds. Successful: %d/%d\n', toc, numel(ok_idx), n_to_run);
    end
end

%% function – moved part of the code here, logic is the same as before
function [oC_i, pC_i, kC_i, RnEijC_i, converged, err_msg] = ...
    run_single_simulation(i, M_const, eps_true_const, rho_true_const, ...
    eps_true_el_const, rho_true_el_const, theta_true_const, knots, ...
    sigma, basis, lnn, kappa_tau, kappa_eps,do_GFT)
    
    oC_i = [];
    pC_i = [];
    kC_i = [];
    RnEijC_i = [];
    converged = false;
    err_msg = [];
    
    rng(123+i)
    try
        Mi = M_const.Value(:,:,i);
        d  = make_data_simulation(Mi, knots, 'cubic');
        % add "true" objects
        d.sigma = sigma;
        d.epsilon_true = eps_true_const;
        d.rho_true     = rho_true_const;
        d.epsilon_true_el = eps_true_el_const.Value(:, i);
        d.rho_true_el     = rho_true_el_const.Value(:, i);
        d.theta_true   = theta_true_const.Value(:, i);
%         d.w_hat_true   = GFT_T_const.Value(:, i);
        
         % knots and n_ij
        kC_i     = d.k;
        RnEijC_i = d.R_nE_ij(:);

        % Estimation
        
        [~, o_i] = GMM_wrapper_gravity(d, kappa_tau, kappa_eps, 'base', basis);
        o_i.title = 'Simulation';

        R1 = 1:size(d.Deriv,2);
        R2 = size(d.Deriv,2)+1:2*size(d.Deriv,2);
        
        % Recover epsilon and rho
        if strcmp(basis,'LogNormal')
            [o_i.epsilon, o_i.rho, o_i.epsilon_el, o_i.rho_el] = recover_elasticities_LN(o_i, d);
        else
            o_i.epsilon_el = d.Deriv * o_i.ests2(R1)';
            o_i.rho_el     = d.Deriv * o_i.ests2(R2)';
            ln       = @(n) evknots_make(d.k,log(n));
            o_i.rho     = exp((ln(exp(d.R_nE_ij))*o_i.est_rho')')';
            o_i.epsilon = exp((ln(exp(d.R_nE_ij))*o_i.est_epsilon')')';

            lnnn                   = @(n) evknots_make(d.k,log(n)) ;
            o_i.rho_est_f         = @(n) exp((lnnn(n)*o_i.est_rho')');                             
            o_i.epsilon_est_f     = @(n) exp((lnnn(n)*o_i.est_epsilon')');

        end
        o_i.theta = (d.sigma-1) * (1 - (1 + o_i.rho_el) ./ o_i.epsilon_el);
        
        % Interpolation using passed lnn
        p_i = function_interp_wrapper_GMM(o_i, d, lnn,do_GFT);
        disp("GFT done")
        
        oC_i = o_i;
        pC_i = p_i;
        converged = true;

    catch ME
        err_msg = sprintf('i=%d\n%s', i, getReport(ME,'extended','hyperlinks','off'));
        converged = false;
    end
end