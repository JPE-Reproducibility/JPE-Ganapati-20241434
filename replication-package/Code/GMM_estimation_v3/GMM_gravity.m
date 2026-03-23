function [ OUT, FE_EST  ] = GMM_gravity( d, epsilons,rhos_bar,kappa_tau,kappa_epsilon,type, W,fix,gamma,gravity )
%% Subroutine for getting out objective function G
% Type designates output type
% type== 0: G*Z
% type== 1: G ( Nx2 vector)
% type== 2: G (2Nx1 vector)
    
    gravity_dim = size(gravity,2);
    d.Z1 = gravity*gamma(1:gravity_dim)'               + d.ltariff_raw ;
    d.Z2 = gravity*gamma(gravity_dim+1:2*gravity_dim)' + d.ltariff_raw ;

    %% Main GMM
    A_EST   = d.Z1                                      - d.X_K * (epsilons') * kappa_epsilon  ;
    B_EST   = d.R_xbar_ij + (d.sigma-1)*kappa_tau.*d.Z2 - d.X_K * (rhos_bar');
    
    % Net Out Estimates
    EST     = [A_EST; B_EST]; 
    FE_EST  = d.FE_PZ_Ki*EST;

    % Recover GMM functions
    G_A     = A_EST - [d.FEA d.FE_C]*FE_EST;
    G_B     = B_EST - [d.FEB d.FE_C]*FE_EST;


    %% Output types - first for GMM; second and third for SE computation
    if type == 0
        % Get Criterion function
        F   = [G_A'*d.ZA_K G_B'*d.ZB_K];
        C   = F*W*F';

        if fix == 1
            C = C + 1000000000*K_ERR;
        end

        OUT = C;
    elseif type == 1
        OUT.G_A = G_A;
        OUT.G_B = G_B;
    elseif type == 2
        OUT = [G_A ; G_B];  
    elseif type == 3

        FEJ_1   = max(d.FEA)'.*FE_EST;  
        FEJ_2   = max(d.FEB)'.*FE_EST;
        f1 = FEJ_1((size(d.R_I,2)*2+1):(size(d.R_I,2)*2+size(d.R_J,2)-1));
        f2 = FEJ_2((size(d.R_I,2)*2+size(d.R_J,2)):(size(FE_EST,1)));
        
        fct_mean = @(x) median(x);
        xi_tilde_epsilon_mean = fct_mean(f1);
        xi_tilde_rho_mean     = fct_mean(f2);
        kappa_r_implied = xi_tilde_epsilon_mean/xi_tilde_rho_mean;
        kappa_tilde_tau     = (1-d.sigma)*kappa_tau;
        kappa_tilde_eps_f   = @(kappa_r) (kappa_tilde_tau*(1+(1-kappa_r)/kappa_r))^-1;
        kappa_tilde_eps_f(kappa_r_implied)

    end
end


