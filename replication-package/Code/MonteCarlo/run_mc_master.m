clear all
close all

%% Master File for Monte Carlo Simulations
% 
% This script handles data generation, GMM estimation, QQ estimation, and plotting
% for different specifications

%% DATA GENERATION ALGORITHM
% 
% 0) Set the functional form(s). Available functional forms are: 
% 'LogNormal', 'MPareto', 'MEstimates'
%
%
% 1.LogNormal – decreasing elasticity
% 2.MPareto – modified Pareto that generates epsilon and G^e as
%  G^e(e) = 1 - (e/e_bar)^(-alpha_e) * (ln(e)/ln(e_bar))^(-gamma_e),
%  gamma_e = 0
% 3.MEstimates - modified Pareto that generates epsilon and G^e as
%  G^e(e) = 1 - (e/e_bar)^(-alpha_e) * (ln(e)/ln(e_bar))^(-gamma_e),
%  gamma_e = 1.5.
%
% 
% 1) For each functional form, we set parameters given the functional form
% using get_fct_form_params.m 
%
% 2) Given the parameters from 1), we use create_elasticity_functions.m
% to generate epsilon, rho, kappa, Q_lnx (sales dist) etc. They are all set
% as function handles
%  
% 3) Given 2) we solve the model using solve_FP.m It always
% uses "parfor" which can be a bit slower than just using "for", but not
% critically. solve_FP.m saves the simulated data files to the specified
% folder
%
%% RUNNING GMM
% 0) Set the specification-estimator map.
%
% 1) Conditional on the data and 0), we run the GMM using by specifying the
% functional form and the estimator to Script_GMM_simulation_ST_alt.m It
% will look for the data file with the functional form in the name
%% RUNNING QQ
% 0) Set the specification-estimator map. Currently we use both 
% n = 1 (n1 option) and n != 1 (nv option) for LogNormal and MLogNormal for testing. 
% Otherwise we set n=1. knots have to always be = 1 and the basis always
% has to be lognormal
%
%
% Note: both Script_GMM_simulation_ST_alt and Script_QQ_simulation_ST_alt
% will now automatically use parfor if more than one core is requested in
% the master file. 
%% RUNNING SPECIFIC SIMULATIONS/ECONOMIES:
% By default, all simulations (1:sims) are run. To run specific simulations:
%
% 1. Single simulation for testing:
%    Script_GMM_simulation_ST_alt(..., 'RunSingle', 5)  % Run only simulation #5
%    Script_QQ_simulation_ST_alt(..., 'RunSingle', 5)
%
% 2. Range of simulations:
%    Script_GMM_simulation_ST_alt(..., 'RunRange', [10 20])  % Run simulations 10-20
%    Script_QQ_simulation_ST_alt(..., 'RunRange', [10 20])
%
% 3. All the simulations (default):
%    Script_GMM_simulation_ST_alt(...)
%    Script_QQ_simulation_ST_alt(...)
% 
% Note: Cannot specify both RunSingle and RunRange simultaneously.
%       For single-core runs (cores=1), progress will be printed sequentially.
%% Plotting
% figures_gen will accept a
% specification-estimator map fo the GMM and a separate
% specification-estimator map for the QQ (because we only plot one QQ)
% % Examples of use:
%
% % If you want individual panels:
% figures_gen(project_path + "output/estimation", ...
%     project_path + "output/plots", "alt", qq_plot, ...
%     'specifications', fct_forms, ...
%     'gmm_estimators', gmm_plot, ...
%     'save_separate', true, ...
%     'save_combined', false);
% 
% % If you want to plot the GFT:
% figures_gen(project_path + "output/estimation", ...
%     project_path + "output/plots", "alt", qq_plot, ...
%     'gmm_estimators', gmm_plot, ...
%     'specifications', fct_forms, ...
%     'plot_gft', true);
%
% % If you want to change the theta (from 0 to 8 is default) limits:
% figures_gen(project_path + "output/estimation", ...
%     project_path + "output/plots", "alt", qq_plot, ...
%     'gmm_estimators', gmm_plot, ...
%     'specifications', fct_forms, ...
%     'ylim_theta', [0, 20]);
%
% % % If you want to change the GFT plot support (from 0 to 8 is default):
% figures_gen(project_path + "output/estimation", ...
%     project_path + "output/plots", "alt", qq_plot, ...
%     'gmm_estimators', gmm_plot, ...
%     'specifications', fct_forms, ...
%     'xlim_gft', [0, 10]);
%
% generate_QQ_plot will accept a simulation id and a pair of countries ids
%
% to plot a simulation 7, countries 3, 8 using estimated and
% true quantiles from MEstimates
% generate_QQ_plot([], output_path + "/data", output_path + "/plots/", ...
%    'MEstimates', 'LogNormal', 1, 'n1', ...
%    [7], {[3,8]});

%% Parameters

project_path = string(fileparts(mfilename('fullpath')));
stack = dbstack('-completenames');
if isempty(stack) || contains(stack(1).file, "LiveEditorEvaluationHelper")
    try
        project_path = string(fileparts(matlab.desktop.editor.getActiveFilename));
    catch
        project_path = "";
    end
end
% command line fallback
if isempty(project_path) || strtrim(project_path) == ""
    project_path = string(pwd);
end

create_simulation = 1;
do_GFT = 1;
run_gmm = 1;
run_qq  = 1;

% Specifications to run
fct_forms = {'LogNormal', 'MEstimates', 'MPareto'};
%fct_forms = {'LogNormal', 'MEstimates', 'MPareto'};

output_path = project_path + "/output";
output_path = "../../Output/";



c = 100;           % countries
sims = 100;        % economies
gamma = 0.1;       % updating parameter in FP
save_graph = 0;
cores = 1;

%% Define specification-estimator maps
% Baseline GMM estimators (for main plots)
gmm_settings_baseline = struct( ...
    'LogNormal', {{struct('basis','Spline','knots',58)}}, ...
    'MPareto', {{struct('basis','Spline','knots',5)}}, ...
    'MEstimates', {{struct('basis','Spline','knots',56)}});

% Additional GMM estimators (for appendix plots)
gmm_settings_additional = struct( ...
    'LogNormal', {{...
        struct('basis','Spline','knots',1), ...
        struct('basis','LogNormal','knots',1)}}, ...
    'MPareto', {{...
        struct('basis','Spline','knots',1), ...
        struct('basis','LogNormal','knots',1)}}, ...
    'MEstimates', {{...
        struct('basis','Spline','knots',1), ...
        struct('basis','LogNormal','knots',1)}});

%% QQ specification-estimator maps
qq_settings = struct( ...
    'LogNormal', {{...
        struct('basis','LogNormal','knots',1, 'force_n','n1')}}, ...
    'MPareto', {{...
        struct('basis','LogNormal','knots',1, 'force_n','n1')}}, ...
    'MEstimates', {{...
        struct('basis','LogNormal','knots',1, 'force_n','n1')}});

%% Plotting specification-estimator maps
% which QQ to plot (pick one per specification)
qq_plot = struct( ...
    'LogNormal', "n1", ...
    'MPareto', "n1", ...
    'MEstimates', "n1");

% GMM estimators for main plots (baseline only)
gmm_plot_main = gmm_settings_baseline;

% GMM estimators for appendix plots (baseline + additional)
gmm_plot_appendix = struct();
for i = 1:length(fieldnames(gmm_settings_baseline))
    specs = fieldnames(gmm_settings_baseline);
    spec = specs{i};
    gmm_plot_appendix.(spec) = [gmm_settings_baseline.(spec), ...
                                  gmm_settings_additional.(spec)];
end

%% Generate the data
if create_simulation == 1
    disp("Generate the data")
    for i = 1:length(fct_forms)
        ff = fct_forms{i};
        
        % Get default parameters for this functional form
        data_params = get_fct_form_params(ff);
        
        % Optional: Override parameters here if needed
        % data_params.sigma_z = 0.7;
        % data_params.mu_f = 0.05;
        
        % Create elasticity functions
        ef = create_elasticity_functions(ff, project_path, output_path, data_params);
        
        % Solve equilibrium
        solve_FP(project_path, output_path, data_params, ef.epsilon, ...
            ef.eps_true, ef.rho, ef.r_true, ef.kappa, ef.H_e, ...
            ef.Q_lnx, ef.p_grid, sims, cores, gamma);
    end
end

%% Run GMM - Baseline
if run_gmm == 1
    disp("Run GMM - Baseline")
    for i = 1:length(fct_forms)
        ff = fct_forms{i};
        
        % Get baseline GMM settings for this specification
        gmm_list = gmm_settings_baseline.(ff);
        
        % Run baseline GMM estimators
        for j = 1:length(gmm_list)
            estimator = gmm_list{j};
            Script_GMM_simulation(project_path, output_path, ff, ...
                {estimator}, cores, sims, do_GFT);
        end
    end
end

%% Run QQ
if run_qq == 1
    disp("QQ Runs")
    
    for i = 1:length(fct_forms)
        ff = fct_forms{i};
        
        % Get QQ settings for this specification
        qq_list = qq_settings.(ff);
        
        % Run QQ for each estimator
        for j = 1:length(qq_list)
            estimator = qq_list{j};
            Script_QQ_simulation(project_path, output_path, ff, ...
                {estimator}, cores, sims, do_GFT)
        end
    end
end

%% Main Plots
for i = 1:length(fct_forms)
    ff = fct_forms{i};
    
    figures_gen(output_path + "/estimation", ...
        output_path + "/plots", "alt", qq_plot, ...
        'specifications', {ff}, ...
        'gmm_estimators', gmm_plot_main.(ff), ...
        'save_separate', true, ...
        'save_combined', false, ...
        'plot_gft', true, ...
        'xlim_gft', [0 3.6], ...
        'ylim_theta', [3, 8], ...
        'gft_tick_delta', 0.3, ...
        'gft_binwidth', 0.03, ...
        'font_size', 15);
end

%% Run GMM - Additional
for i = 1:length(fct_forms)
    ff = fct_forms{i};
    
    % Get additional GMM settings for this specification
    gmm_list = gmm_settings_additional.(ff);
    
    % Run additional GMM estimators
    for j = 1:length(gmm_list)
        estimator = gmm_list{j};
        Script_GMM_simulation(project_path, output_path, ff, ...
            {estimator}, cores, sims, 0);
    end
end

%% Appendix Plots (with all GMM estimators)
for i = 1:length(fct_forms)
    ff = fct_forms{i};
    
    figures_gen(output_path + "/estimation", ...
        output_path + "/plots/additional specs", "alt", qq_plot, ...
        'specifications', {ff}, ...
        'gmm_estimators', gmm_plot_appendix.(ff), ...
        'save_separate', true, ...
        'save_combined', false, ...
        'ylim_theta', [3 30],...
        'file_suffix', 'apx', ...
        'legend_position', 'northeast', ...
        'font_size', 15);
end


%% QQ plot
generate_QQ_plot([], output_path + "/data", output_path + "/plots/", ...
    'MEstimates', 'LogNormal', 1, 'n1', ...
    [1], {[1,4]});
