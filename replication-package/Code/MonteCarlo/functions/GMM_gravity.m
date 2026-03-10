function [ OUT, FE_EST  ] = GMM_gravity( d, epsilons,rhos,kappa_tau,kappa_epsilon,type, W )
%% Subroutine for getting out objective function G
% Type designates output type
% type== 0: G*Z
% type== 1: G ( Nx2 vector)
% type== 2: G (2Nx1 vector)
    d.Z = d.ltariff_raw ;


    %% Main GMM
    A_EST   = d.Z                                      - d.X_K * (epsilons') * kappa_epsilon  ;
    B_EST   = d.R_xbar_ij + (d.sigma-1)*kappa_tau.*d.Z - d.X_K * (rhos');


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
        OUT = C;
    elseif type == 1
        OUT.G_A = G_A;
        OUT.G_B = G_B;
    elseif type == 2
        OUT = [G_A ; G_B];  
    end

end


