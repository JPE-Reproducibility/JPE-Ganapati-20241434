function l = regQQ(d, Q_lnx_all, p_grid, basis, knots, force_n)
    % Modified to always store data needed for QQ plotting
    
    %Implement an FE regression and obtain rho and epsilon
    i_idx          = d.I;
    j_idx          = d.J;
    n_ij           = exp(d.R_nE_ij);
    c              = sqrt(size(i_idx, 1));
    pp = numel(p_grid);
    
    %Set up FEs
    FE_Q = sub2ind([c c], i_idx, j_idx);
    FE_ij = sparse(1:numel(FE_Q), FE_Q, 1, numel(FE_Q), c^2);
    FE_ij = repmat(FE_ij, pp, 1);
    G = FE_ij;
    cnt = full(FE_ij' * ones(size(Q_lnx_all,1),1));
    
    %Group means for y
    Ybar      = (G' * Q_lnx_all) ./ cnt;          % y means
    y_within  = Q_lnx_all - G * Ybar;             % y
    
    if strcmp(basis,'LogNormal')
        %group means for x
        if strcmp(force_n, "n1")
            n_ij_alt  = ones(size(n_ij));
            m_ij      = norminv(1 - n_ij_alt*((1 - p_grid))); %testing
        else
            m_ij      = norminv(1 - n_ij*((1 - p_grid))); %testing
        end
        m_ij      = m_ij(:);
        Mbar      = (G' * m_ij) ./ cnt;
        m_within  = m_ij - G * Mbar;
        
        % run a demeaned reg w/ const
        X1 = [ones(c^2*pp, 1), m_within];
        [b, ~, ~, ~, stats] = regress(y_within, X1);
        l.R2_within  = stats(1);   % K x 1
        l.beta = b(2);
        l.beta0 = b(1);
        o.est_sigma = l.beta;
        
        [l.epsilon, l.rho, l.epsilon_el, l.rho_el] = recover_elasticities_LN(o, d);
        
    elseif strcmp(basis, 'Spline') && knots == 1
        %group means for x
        m_ij      = log(1 - p_grid) + 0*log(n_ij); %testing
        Mbar      = (G' * m_ij) ./ cnt;                  % m means
        m_within  = m_ij - G * Mbar;                     % m
        
        % run a demeaned reg w/ const
        X1 = [ones(c^2*pp, 1), m_within];
        [b, ~, ~, ~, stats] = regress(y_within, X1);
        l.R2_within  = stats(1);   % K x 1
        l.beta = b(2);
        l.beta0 = b(1);
        l.rho = l.beta;
        l.epsilon = l.beta;
    end
    
    % Always store the full demeaned data for plotting
    % Users can extract specific i-j pairs later
    l.y = Q_lnx_all;
    l.X = [ones(c^2*pp, 1), m_ij];
    l.y_within = y_within;  % Full demeaned Y data [c^2*pp x 1]
    l.X1 = X1;              % Full demeaned X data [c^2*pp x 2]
    l.c = c;                % Number of countries
    l.pp = pp;              % Number of quantiles
end