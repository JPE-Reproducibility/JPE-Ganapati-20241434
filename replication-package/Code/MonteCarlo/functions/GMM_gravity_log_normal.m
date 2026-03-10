function [ OUT, FE_EST  ] = GMM_gravity_log_normal( d, ~, n_sigma, kappa_tau,kappa_epsilon,type, W )
%% Subroutine for getting out objective function G
% Type designates output type
% type== 0: G*Z
% type== 1: G ( Nx2 vector)
% type== 2: G (2Nx1 vector)
    
    d.Z = d.ltariff_raw ;

%     % Log Normal Functions     
%     lne = @(lnn) n_sigma * norminv(1-exp(lnn), 0, n_sigma);
%     r = @(n) 1./n .* arrayfun(@(l,u) integral( @(n) exp(n_sigma * norminv(1-n, 0, n_sigma)) ,l,u),zeros(size(n)),n);
%     lnr = @(lnn) log(r(exp(lnn)));

    do_interpolation = 1 ;
    if do_interpolation == 1
        % Interpolation Function
        lower = min(min(d.X_K)-1,-4);
        u_vals = exp([lower:.01:-.0001 -.00001]);
        integrals = zeros(size(u_vals)); 
        d.X_K(d.X_K>log(max(u_vals))) = log(max(u_vals));
        
        for i = 1:length(u_vals)
            integrals(i) = integral(@(n) exp(n_sigma * norminv(1 - n, 0, 1)), 0, u_vals(i)); %,'AbsTol',1e-20,'RelTol',1e-10
        end

        r = @(n) 1./n .*  interp1(u_vals, integrals, n, 'linear');
    else
        r = @(n) 1./n .* arrayfun(@(l,u) integral( @(n) exp(n_sigma * norminv(1-n, 0, 1)) ,l,u),zeros(size(n)),n);
    end

    % This is due to the lack of precision in the numerical integration
    d.X_K(d.X_K<-23 ) = -23;

    lne = @(lnn) n_sigma * norminv(1-exp(lnn), 0, 1);
    lnr = @(lnn) log(r(exp(lnn)));

    %% Main GMM
    A_EST   = d.Z                                      - lne(d.X_K) * kappa_epsilon  ;
    B_EST   = d.R_xbar_ij + (d.sigma-1)*kappa_tau.*d.Z - lnr(d.X_K);

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


