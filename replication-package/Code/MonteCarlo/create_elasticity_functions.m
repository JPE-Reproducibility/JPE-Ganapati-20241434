function out = create_elasticity_functions(fct_form, project_path, output_path, params)

arguments
    fct_form string {mustBeMember(fct_form, ["LogNormal", ...
        "MPareto", "MEstimates"]),  mustBeNonempty}
    project_path string {mustBeNonempty}
    output_path string {mustBeNonempty}
    params struct
end
%paths
data_path = output_path + '/data';
addpath(data_path, '-frozen')
addpath(project_path + '/functions', '-frozen')
%params
sigma = params.sigma;
nu_ij = params.nu_ij;
nu    = params.nu; % used for Mlognormal only
eps   = 1e-3;
% sigma is the elasticity of substitution parameter
% nu is the dispersion parameter in the sales distribution
% eps is the step size in the elasticity computation

%define grid for computation
% log grid for n_ij
ln_grid = linspace(-20, -.0001, 2000)'; 
% percentiles of sales distribution to be computed
p_grid = (1:99)/100;                     % 1..99 percent

     if strcmp(fct_form,'LogNormal') == 1
         % Lognormal
        
        % extensive margin elasticity function
        z_f         = @(n) exp(sqrt(2)*nu_ij*erfinv(2*(1-n)-1));
        %z_f        = @(n) interp1(n_grid, z_f_grid, n, 'pchip', 'extrap');
        e_f         = @(n) (n);
        epsilon     = @(ln) e_f(z_f(exp(ln)));

        %elasticity
        lne         = @(lnn) log(epsilon(lnn));
        eps_true    = @(n) (lne(log(n)+eps)-lne(log(n)))./(eps);
        
        r_int_f     = @(n) integral(@(n2) e_f(z_f(n2)),0,n);
        r_f         = @(n) 1./n .* arrayfun(@(l,u) integral( @(n2) e_f(z_f(n2)) ,l,u),zeros(size(n)),n);
        r_ex        = @(n) r_f(n);
        r_f_grid    = arrayfun(@(ln) r_f(exp(ln)), ln_grid);
        
        % intensive margin elasticity function
        rho         = @(n) interp1(ln_grid, r_f_grid, log(n), 'spline', 'extrap'); % why are we interpolating. i need a good reason.
        % rho       = @(n) 1./n .* arrayfun(@(l,u) integral( @(n2) e_f(z_f(n2)) ,l,u),zeros(size(n)),n);
        % rho       = @(n) 1./n .* arrayfun(@(l,u) integral( @(n2) e_f(z_f(n2)) ,l,u),zeros(size(n)),n);

        rho_bar     = @(n) e_f(z_f(n));

        % elasticity
        lnr         = @(lnn) log(rho(exp(lnn)));
        r_true      = @(n) (lnr(log(n)+eps)-lnr(log(n)))./(eps);
        
        %kappa
        ratio       = @(n) rho_bar(n)./epsilon(log(n));
        kappa       = @(n) arrayfun(@(lb,ub) integral(@(k) ratio(k),lb,ub), zeros(size(n)), n);
        kappa_grid  = arrayfun(@(ln) kappa(exp(ln)), ln_grid);
        kappa_int   = @(n) interp1(ln_grid, kappa_grid, log(n), 'spline', 'extrap');
        % H^e
        
        % I MADE A CHANGE HERE from logncdf directly
        H_e         = @(lne) logncdf(exp(lne), 0, nu_ij);
        

        %Sales distribution
        m_ij        = @(n) norminv(1 - ((1 - p_grid))); 
        %Q1_grid     = m_ij(exp(ln_grid)); % double grid
        %Q1_int      = @(n) interp1(ln_grid, Q1_grid, log(n), 'spline', 'extrap');
        Q_lnx       = @(R_bar, n) log(R_bar) + nu_ij * m_ij(n); %nu_ij defined above

    elseif strcmp(fct_form,'MPareto') == 1
        
        
        e_bar       = exp(1);
        alpha_e     = 5./(sigma-1) - 1;
        gamma_e     = 0;
        
        %% Functions
    
        k1 = -(sigma - 1)/2000;
        k2 = (2000 - sigma + 1)/2000;
        
    
        %e_power   = k1;
        rho_power = k1;
        n_power   = k2;
         
        % get true theta, force it to 5
        epsilon     = @(ln) epsilon_f(exp(ln), e_bar, alpha_e, gamma_e);
        lne         = @(lnn) log(epsilon_f(exp(lnn), e_bar, alpha_e, gamma_e));
        eps_true    = @(n) (lne(log(n)+eps)-lne(log(n)))./(eps);

        
        % % intensive margin elasticity function, set rho_bar and rho close to 1
        rho_bar     = @(n) (rho_power+1)*n.^rho_power./n_power;
        rho         = @(n) n.^(rho_power)./n_power;

        lnr         = @(lnn) log(rho(exp(lnn)));
        r_true      = @(n) (lnr(log(n)+eps)-lnr(log(n)))./(eps);
        
        % kappa
        ratio       = @(ln) (rho_bar(exp(ln))./epsilon(ln));
        kappa       = @(ln) [0; cumsum(exp(ln(1:end-1)).*ratio(ln(1:end-1)) .* diff(ln))];
        kappa_grid  = kappa(ln_grid);
        kappa_int   = @(n) interp1(ln_grid, kappa_grid, log(n), 'spline', 'extrap');
        
        % H^e
        H_e         = @(lne) generate_H(exp(lne), e_bar, alpha_e, gamma_e);
        %inv_H_e = @(n) interp1(1 - n_sorted, e_sorted, n, 'linear', 'extrap');
        
        rho_var     = @(n) r_true(n);
        
        % Sales distribution
        p_col   = p_grid(:);                       % ppx1
        p_inv   = norminv(p_col).';                % 1xpp
        x_grid  = @(n, rb) log(rho(n)) + log(1+rho_var(n)) + log(rb) + nu_ij.*p_inv;
    
        lnH  = @(xx, s, rb) normcdf( ...
            ( xx - log(rho(s)) - log(1+rho_var(s)) - log(rb)) ./ nu_ij );
    
        Q_lnx  = @(Rbar, n) inv_unique_interp(lnH, x_grid, p_grid, n, Rbar);

    elseif strcmp(fct_form,'MEstimates') == 1
        
        
        e_bar       = exp(1);
        alpha_e     = 5./(sigma-1) - 1;
        gamma_e     = 1.5;
        
        %% Functions
    
        k1 = -(sigma - 1)/2000;
        k2 = (2000 - sigma + 1)/2000;
        
    
        %e_power   = k1;
        rho_power = k1;
        n_power   = k2;
         
        % get true theta, force it to 5
        epsilon     = @(ln) epsilon_f(exp(ln), e_bar, alpha_e, gamma_e);
        lne         = @(lnn) log(epsilon_f(exp(lnn), e_bar, alpha_e, gamma_e));
        eps_true    = @(n) (lne(log(n)+eps)-lne(log(n)))./(eps);

        
        % % intensive margin elasticity function, set rho_bar and rho close to 1
        rho_bar     = @(n) (rho_power+1)*n.^rho_power./n_power;
        rho         = @(n) n.^(rho_power)./n_power;

        lnr         = @(lnn) log(rho(exp(lnn)));
        r_true      = @(n) (lnr(log(n)+eps)-lnr(log(n)))./(eps);
        
        % kappa
        ratio       = @(ln) (rho_bar(exp(ln))./epsilon(ln));
        kappa       = @(ln) [0; cumsum(exp(ln(1:end-1)).*ratio(ln(1:end-1)) .* diff(ln))];
        kappa_grid  = kappa(ln_grid);
        kappa_int   = @(n) interp1(ln_grid, kappa_grid, log(n), 'spline', 'extrap');

        % H^e
        H_e         = @(lne) generate_H(exp(lne), e_bar, alpha_e, gamma_e);
        %inv_H_e = @(n) interp1(1 - n_sorted, e_sorted, n, 'linear', 'extrap');
        
        rho_var     = @(n) r_true(n);
        
        % Sales distribution
        p_col   = p_grid(:);                       % ppx1
        p_inv   = norminv(p_col).';                % 1xpp
        x_grid  = @(n, rb) log(rho(n)) + log(1+rho_var(n)) + log(rb) + nu_ij.*p_inv;
    
        lnH  = @(xx, s, rb) normcdf( ...
            ( xx - log(rho(s)) - log(1+rho_var(s)) - log(rb)) ./ nu_ij );
    
        Q_lnx  = @(Rbar, n) inv_unique_interp(lnH, x_grid, p_grid, n, Rbar);
    end

out.rho = rho;
out.r_true = r_true;
out.H_e = H_e;
out.epsilon = epsilon;
out.eps_true = eps_true;
out.kappa = kappa_int;
out.Q_lnx = Q_lnx;
out.p_grid = p_grid;

end