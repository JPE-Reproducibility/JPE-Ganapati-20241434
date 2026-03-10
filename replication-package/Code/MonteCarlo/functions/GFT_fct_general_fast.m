function [W_hat_GFT] = GFT_fct_general_fast(epsilon,rho,sigma,X_ij,n_ij,trigger)

    % For Testing
    %     n_mu          = o.ests2(1);
    %     n_sigma       = o.ests2(2); 
    %     z_LN          = @(n) exp(n_mu + sqrt(2)*n_sigma*erfinv(2*(1-n)-1));
    %     e_LN          = @(n) (n);
    %     epsilon       = @(n) e_LN(z_LN(n));
    %     rho           = @(n) e_LN(z_LN(n));

    C                = size(n_ij,1);
    x_ij             = X_ij ./ repmat(sum(X_ij,1),C,1);
    rho_bar          = @(n) 1./n .* arrayfun(@(l,u) integral( @(n2) rho((n2)) ,l,u),zeros(size(n)),n);
    ratio            = @(n) rho(n)./epsilon(n);


    do_pre_interp = 0;
    do_interpolation = 0;
    if nargin <= 6
        trigger = 0;
        buffer = 0;
    else
        if trigger == 1
            do_pre_interp = 0;
            buffer = 0;
        end
    end
    

    if do_pre_interp == 1
        disp("Pre-Calc Running")
        % One-time setup: Create interpolation table
        n_grid = logspace(log(min(min(n_ij))), 0, 500); % Log-spaced from 0.0001 to 1
        
        % Pre-compute integrals
        top_interp = zeros(size(n_grid));
        bot_interp = zeros(size(n_grid));
    
        % First interval
        top_interp(1) = integral(ratio, 0, n_grid(1));
        bot_interp(1) = integral(rho, 0, n_grid(1));
        
        % Cumulative sum: only integrate small segments
        for i = 2:length(n_grid)
            top_interp(i) = top_interp(i-1) + integral(ratio, n_grid(i-1), n_grid(i));
            bot_interp(i) = bot_interp(i-1) + integral(rho, n_grid(i-1), n_grid(i));
        end
    
        disp("Interpolation Setup")
        % Create interpolants (griddedInterpolant is faster than interp1)
        top_interpolant = griddedInterpolant(n_grid, top_interp, 'pchip');
        bot_interpolant = griddedInterpolant(n_grid, bot_interp, 'pchip');
        %     eps_interpolant = griddedInterpolant(n_grid, eps_interp, 'pchip');
    
        disp("Interpolation Running")
        % Fast evaluation (vectorized)
        top1 = top_interpolant(n_ij);
        bot1 = bot_interpolant(n_ij) ./ arrayfun(epsilon,n_ij);

    elseif do_pre_interp == 2

        disp("Interpolation Setup")
        % Extract and sort unique values from n_ij
        n_unique = unique(n_ij(:));
        n_unique = [0; n_unique(:)]; % Add 0 as starting point
        
        % Compute cumulative integrals only for unique values
        top_cumulative = zeros(size(n_unique));
        bot_cumulative = zeros(size(n_unique));
        
        for i = 2:length(n_unique)
            disp(i)
            top_cumulative(i) = top_cumulative(i-1) + integral(ratio, n_unique(i-1), n_unique(i));
            bot_cumulative(i) = bot_cumulative(i-1) + integral(rho, n_unique(i-1), n_unique(i));
        end

        disp("Interpolation Running")
        % Map back to original R×C structure using interpolation
        top1 = interp1(n_unique, top_cumulative, n_ij, 'pchip');
        bot1 = interp1(n_unique, bot_cumulative, n_ij, 'pchip') ./ epsilon(n_ij);
    else
        disp("Pre-Calc Running - Integrals")
        top1             = arrayfun(@(lb,ub) integral(ratio,lb,ub),zeros(C),n_ij);
        bot1             = arrayfun(@(lb,ub) integral(rho,lb,ub)  ,zeros(C),n_ij)./arrayfun(epsilon,n_ij);
    end


    pre_integrals    = top1./bot1; 
    intial_guess     = ones(size(n_ij,1),1)+.01;

    if do_interpolation == 1
        % Interpolation Function
        lower =-10;
        u_vals = exp([lower:.01:-.0001 -.00001]);
        top_integrals = zeros(size(u_vals)); 
        bot_integrals = zeros(size(u_vals)); 

        for i = 1:length(u_vals)
            top_integrals(i) = integral(ratio, 0, u_vals(i)); %,'AbsTol',1e-20,'RelTol',1e-10
            bot_integrals(i) = integral(rho, 0, u_vals(i)); %,'AbsTol',1e-20,'RelTol',1e-10
        end

        top_fint = @(lb,n) interp1(u_vals, top_integrals, n, 'linear');
        bot_fint = @(lb,n) interp1(u_vals, bot_integrals, n, 'linear');

    else
    %         % this is the slowest step - approximations can help!
    %         % closed forms for pareto
        top_fint = @(lb,ub) integral(ratio,lb,ub); 
        bot_fint = @(lb,ub) integral(rho,lb,ub);   
    end
 
    FGFT             = @(n_ii_hat) GFT_pre(n_ii_hat,X_ij,epsilon,rho,rho_bar,ratio,n_ij,x_ij,pre_integrals,top_fint,bot_fint);
%     [n_ii_hat]       = fsolve(FGFT,intial_guess,optimoptions('fsolve','Display','iter'));

    disp("GFT Running")
    n_ii_hat = ones(size(intial_guess));
    for i=1:C
        fprintf('%3.0f,',i)
        FGFT1             = @(n_ii_hat) GFT_pre1(i,n_ii_hat,X_ij,epsilon,rho,rho_bar,ratio,n_ij,x_ij,pre_integrals,top_fint,bot_fint);
        [X,~,EXITFLAG] ...
                          = fsolve(FGFT1,intial_guess(i),optimoptions('fsolve','Display','none', 'OptimalityTolerance', 1e-11, 'FunctionTolerance', 1e-11, 'MaxFunctionEvaluations',10000, 'MaxIterations',10000));
        n_ii_hat(i)       = X;  
        
        %{
        % This is about crazy outliers.
        if EXITFLAG < 0
            W_hat_GFT = NaN;
            return 
        end
        %}
    end


    disp("GFT Aggregated")

    [~,N_hat_LN]     = FGFT(n_ii_hat);
    elasticity       = -1/(sigma-1);
    W_hat_GFTLN_1    = diag(x_ij).^elasticity;
    W_hat_GFTLN_2    = (n_ii_hat.*N_hat_LN).^elasticity;
    W_hat_GFTLN_3    = ( arrayfun(rho_bar,n_ii_hat.*diag(n_ij))./arrayfun(rho_bar,diag(n_ij))    ).^elasticity;
    W_hat_GFT        = W_hat_GFTLN_1 .* W_hat_GFTLN_2 .* W_hat_GFTLN_3;

end

