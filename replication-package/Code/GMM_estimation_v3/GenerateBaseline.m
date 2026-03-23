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
    save 'baseline_estimate' 'o'


end