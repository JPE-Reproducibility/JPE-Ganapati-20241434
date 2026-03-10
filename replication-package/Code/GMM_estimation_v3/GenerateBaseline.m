function  GenerateBaseline
% Recover Baseline Estimate

%% Using Estimates
    addpath('../Estimation/GMM_estimation_v3')
    data_dir = '../../Data/Int/WIOD_sampleB/';

%% Using Estimates - knot spline
    basefile = '2012_Tetiuw.csv'; name = '2012_Teti';  grav_spec = 'all';
    load('gamma_gravity.mat');
    % gamma_gravity =    [0.3520   -0.0214   -0.0189   -0.2324   -0.1502    0.2812   -0.0834   -0.0015   -0.0576   -0.0471];
    kappa_tau = 1;
    kappa_f   = 0;
    sigma     = 3.2;
    kappa_epsilon = 1/((sigma-1)*kappa_tau+kappa_f);
    knots = 3;
    M_base = csvread(strcat(data_dir,basefile));
    d = make_data_combo(M_base,knots,'cubic','','',0,grav_spec);    d.sigma = sigma;  
    [~,o] =  GMM_wrapper_gravity( d, kappa_tau, kappa_epsilon,'base',0,gamma_gravity,grav_spec );
    o.k = d.k;
    %save '../Code/Parametric Algorithm/baseline_estimate' 'o'
    save 'baseline_estimate' 'o'

    %{
    lnn     = @(n) evknots_make(o.k,log(n)) ;
    lnd     = @(n) evknots_make_deriv(o.k,log(n)) ;
    rho     = @(n) exp((lnn(n)*o.est_rho')');                             
    rho_bar = @(n) exp((lnn(n)*o.est_rho')') .* ((lnd(n)*o.est_rho')'+1); 
    epsilon = @(n) exp((lnn(n)*o.est_epsilon')');
    ratio   = @(n) rho_bar(n)./epsilon(n);
    %}

%     C       = size(n_ij,1);
%     top1 = arrayfun(@(lb,ub) integral(ratio,lb,ub),zeros(C),n_ij);
%     bot1 = arrayfun(@(lb,ub) integral(rho,lb,ub)  ,zeros(C),n_ij)./arrayfun(epsilon_bar,n_ij);
%     sum(top1./bot1,[1,2])
% 
%     precalculated_integrals = top1./bot1; % Also Equal to (rho_power+1)/(rho_power+1-e_power);

end