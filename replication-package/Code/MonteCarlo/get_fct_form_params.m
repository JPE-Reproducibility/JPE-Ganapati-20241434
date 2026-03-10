function params = get_fct_form_params(fct_form)
    % Get default data generation parameters for a given functional form
    % 
    % Usage:
    %   params = get_data_params('Pareto', 3.2, 1.2)
    
        arguments
            fct_form string {mustBeMember(fct_form, ["LogNormal", ...
        "MPareto", "MEstimates"]),  mustBeNonempty}
    end
    
    % Common parameters
    params = struct();
    params.fct_form = fct_form;
    params.c     = 100;       % number of countries
    params.sigma = 3.2;       % elasticity of substitution
    params.nu_ij = 1.2;       % dispersion parameter in sales distribution
    params.nu = 1.8;       % dispersion parameter in lognormal
    
    % Set defaults based on functional form
    switch fct_form
            
        case 'LogNormal'
            % Log-Normal (Decreasing Elasticity)
            params.mu_f = 0.3;
            params.sigma_z = 0.35;
            params.sigma_r = 0.000;
            params.sigma_f = 0.000;
            params.kappa_r = -2.2;
            params.kappa_f = 0;
            params.kappa_r_home = 0;
            params.kappa_f_home = 0;
            params.mu_z = 0;
            params.mu_r = -6.6;

        case 'MPareto'
            % Constant Elasticity
            params.mu_f = 5.5;
            params.sigma_z = 0.6;
            params.sigma_r = 0.0;
            params.sigma_f = 0;
            params.kappa_r = -2.2;
            params.kappa_f = 0;
            params.kappa_r_home = 0.0;
            params.kappa_f_home = 0.0;
            params.mu_r = -11;
            params.mu_z = 0;
        case 'MEstimates'
            % Increasing Elasticity

            params.mu_f = 9;
            params.sigma_z = 0.45;
            params.sigma_r = 0;
            params.sigma_f = 0;
            params.kappa_r = -2.2;
            params.kappa_f = 0;
            params.kappa_r_home = 0.0;
            params.kappa_f_home = 0.0;
            params.mu_r = -0.6;
            params.mu_z = 0;
    end
end