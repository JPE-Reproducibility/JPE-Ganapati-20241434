## Filepaths Analysis Details

**/Users/florianoswald/actions-runner/_work/JPE-Ganapati-20241434/JPE-Ganapati-20241434/replication-package/Code/GMM_estimation_v3/Master_script_gravity_v3.m**

- Line 4, unix : data_dir = '../../Data/Int/WIOD_sampleB/';
- Line 34, unix : p_overlay_elasticity( d,o,'theta'       ,'../../Output/full_SE' );
- Line 49, unix : p_overlay_elasticity( d,o,'epsilon_CC'   ,'../../Output/F2a_base_' );
- Line 50, unix : p_overlay_elasticity( d,o,'rho_CC'       ,'../../Output/F2a_base_' );
- Line 51, unix : p_overlay_elasticity( d,o,'extensive_CC'       ,'../../Output/F2a_base_' );
- Line 52, unix : h01 = p_overlay_elasticity( d,o, 'theta', '../../Output/F2a_base_'   );
- Line 60, unix : p_overlay_elasticity( dIV,oIV,'theta'       ,'../../Output/F5_IV_' );
- Line 95, unix : p_overlay_elasticity( d2q_Wo,o2q_Wo,'theta_split'       ,'../../Output/test' );
- Line 102, unix : p_overlay_elasticity( dIVs,oIVs,'theta_split'       ,'../../Output/F5_IVo_' );
- Line 187, unix : h00 = p_overlay_elasticity( p_lin,o_lin, 'theta', '../../Output/test'   );
- Line 208, unix : h03 = p_overlay_elasticity( d4q_Wod,o4q_Wod, 'theta', '../../Output/test'   );
- Line 210, unix : p_overlay_elasticity( d4q_Wod,o4q_Wod,  'theta_quad_nose'  ,'../../Output/F3d1_Cross_' );
- Line 211, unix : p_overlay_elasticity( d4q_Wod,o4q_Wod,  'theta_quadA'  ,'../../Output/F3d1_Cross_' );
- Line 212, unix : p_overlay_elasticity( d4q_Wod,o4q_Wod,  'theta_quadB'  ,'../../Output/F3d1_Cross_' );
- Line 213, unix : p_overlay_elasticity( d4q_Wod,o4q_Wod,  'extensive_quad'  ,'../../Output/F3d1_Cross_' );
- Line 223, unix : h02 = p_overlay_elasticity( d4q_Wod_l,o4q_Wod_l, 'theta_quad', '../../Output/test'   );
- Line 250, unix : saveas(gcf,'../../Output/theta_dist_all2','epsc')
- Line 312, unix : obasekf5.title = '\kappa_e = 1/2';

**/Users/florianoswald/actions-runner/_work/JPE-Ganapati-20241434/JPE-Ganapati-20241434/replication-package/Code/ColombiaTests/QuantilesColombia.m**

- Line 44, unix : saveas(gcf,'../../Output/QQ_Colombia_Baseline','epsc')
- Line 46, unix : sgtitle({'Baseline Economy fitted to Colombia-USA Trade',' n_{ij} = exporters to USA/total Exporters'})
- Line 47, unix : saveas(gcf,'../../Output/QQ_Colombia_Baseline','pdf')
- Line 81, unix : saveas(gcf,'../../Output/QQ_Colombia_LogNormal','pdf')
- Line 114, unix : saveas(gcf,'../../Output/QQ_Colombia_LogNormal_Selection','pdf')
- Line 145, unix : saveas(gcf,'../../Output/QQ_Colombia_Pareto','pdf')
- Line 183, unix : saveas(gcf,'../../Output/QQ_Colombia_Quantiles','pdf')
- Line 195, unix : saveas(gcf,'../../Output/QQ_Colombia_Quantiles_BaselineP','epsc')
- Line 205, unix : saveas(gcf,'../../Output/QQ_Colombia_Quantiles_LNP','epsc')

**/Users/florianoswald/actions-runner/_work/JPE-Ganapati-20241434/JPE-Ganapati-20241434/replication-package/Code/GMM_estimation_v3/CC_se_thetae.m**

- Line 38, unix : grad(i) = (1/Deriv_eps(i).^2);

**/Users/florianoswald/actions-runner/_work/JPE-Ganapati-20241434/JPE-Ganapati-20241434/replication-package/Code/GMM_estimation_v3/GMM_wrapper_gravity.m**

- Line 65, unix : Lambda  = Lambda/d.N;
- Line 128, unix : Lambda  = Lambda/d.N;

**/Users/florianoswald/actions-runner/_work/JPE-Ganapati-20241434/JPE-Ganapati-20241434/replication-package/Code/GMM_estimation_v3/hs_routine_gravity.m**

- Line 55, unix : p_overlay_elasticity( d_hs,o_hs,'theta'       ,['../../Output/F6_HS_' int2str(hs)] ,strcat(title))
- Line 56, unix : p_overlay_elasticity( d_hs,o_hs,'epsilon'       ,['../../Output/F6_HS_' int2str(hs)] ,strcat('Extensive: ',title))
- Line 57, unix : p_overlay_elasticity( d_hs,o_hs,'rho'       ,['../../Output/F6_HS_' int2str(hs)] ,strcat(title))
- Line 58, unix : p_overlay_elasticity( d_hs,o_hs,'extensive'       ,['../../Output/F6_HS_' int2str(hs)] ,strcat(title))
- Line 66, unix : p_overlay_elasticity( d_split_hs,o_split_hs,'theta_split'       ,['../../Output/F6_HS_' int2str(hs)] ,strcat('\theta: ',title))

**/Users/florianoswald/actions-runner/_work/JPE-Ganapati-20241434/JPE-Ganapati-20241434/replication-package/Code/LogCorrected/generate_pdf_entry.m**

- Line 61, windows : fprintf('Integral of PDF (γ^e = %d) over [%.2f, 10]: %.4f\n', ...

**/Users/florianoswald/actions-runner/_work/JPE-Ganapati-20241434/JPE-Ganapati-20241434/replication-package/Code/MonteCarlo/run_mc_master.m**

- Line 52, unix : %% RUNNING SPECIFIC SIMULATIONS/ECONOMIES:

**/Users/florianoswald/actions-runner/_work/JPE-Ganapati-20241434/JPE-Ganapati-20241434/replication-package/Code/LogCorrected/plot_theta_e_vs_log_n.m**

- Line 11, unix : % Formula: theta^e(n) = alpha^e + gamma^e/ln(epsilon(n))

**/Users/florianoswald/actions-runner/_work/JPE-Ganapati-20241434/JPE-Ganapati-20241434/replication-package/Code/p05_Sample_Creation_HS2.do**

- Line 107, unix : gen tariff_use = simpleAHS_w/100
- Line 111, unix : gen tariff_useUW = simpleAHS_uw/100
- Line 115, unix : gen Ttariff_useUW = TsimpleAHS_uw/100

**/Users/florianoswald/actions-runner/_work/JPE-Ganapati-20241434/JPE-Ganapati-20241434/replication-package/Code/MonteCarlo/functions/GMM_wrapper_gravity.m**

- Line 42, unix : Lambda  = Lambda/d.N;
- Line 102, unix : Lambda  = Lambda/d.N;

**/Users/florianoswald/actions-runner/_work/JPE-Ganapati-20241434/JPE-Ganapati-20241434/replication-package/Code/LogCorrected/plot_theta_combined.m**

- Line 17, unix : %   theta^e(n) = alpha^e + gamma^e/ln(epsilon(n))

**/Users/florianoswald/actions-runner/_work/JPE-Ganapati-20241434/JPE-Ganapati-20241434/replication-package/Code/GMM_estimation_v3/Simulation_Data_v10.m**

- Line 3, unix : data_dir = '../../Data/Int/WIOD_sampleB/';
- Line 31, unix : z              = unif.^(-1/theta);
- Line 139, unix : z           = unif_EKK.^(-1/theta_EKK);
- Line 151, unix : elast_r_EKK =     (r_vec/iter)./cumtrapz(unif_EKK,(r_vec/iter))./unif_EKK-1;
- Line 152, unix : e               = e_vec/iter;
- Line 187, unix : saveas(gcf,  '../../Output/Epsilon_v8','epsc')
- Line 206, unix : saveas(gcf,  '../../Output/Rho_v8','epsc')
- Line 224, unix : saveas(gcf,  '../../Output/Theta_v8','epsc')
- Line 245, unix : saveas(gcf,  '../../Output/Extensive_v8','epsc')
- Line 264, unix : saveas(gcf,  '../../Output/Intensive_v8','epsc')
- Line 283, unix : saveas(gcf,  '../../Output/Thetac_v8','epsc')
- Line 356, unix : saveas(gcf,  '../../Output/Theta_overlay_v8','pdf')
- Line 357, unix : saveas(gcf,  '../../Output/Theta_overlay_v8','epsc')
- Line 376, unix : saveas(gcf,  '../../Output/Theta_overlay2_v8','pdf')
- Line 377, unix : saveas(gcf,  '../../Output/Theta_overlay2_v8','epsc')
- Line 397, unix : saveas(gcf,  '../../Output/Rho_overlay_v8','epsc')
- Line 398, unix : saveas(gcf,  '../../Output/Rho_overlay_v8','pdf')
- Line 419, unix : saveas(gcf,  '../../Output/Epsilon_overlay_v8','epsc')
- Line 420, unix : saveas(gcf,  '../../Output/Epsilon_overlay_v8','pdf')
- Line 440, unix : saveas(gcf,'../../Output/Extensive_overlay_v8','epsc')
- Line 441, unix : saveas(gcf, '../../Output/Extensive_overlay_v8','pdf')
- Line 467, unix : saveas(gcf,'../../Output/Intensive_overlay_v8','epsc')
- Line 468, unix : saveas(gcf,'../../Output/Intensive_overlay_v8','pdf')

**/Users/florianoswald/actions-runner/_work/JPE-Ganapati-20241434/JPE-Ganapati-20241434/replication-package/Code/LinearCF_v3/LinearCF_wrapper_v3.m**

- Line 22, unix : weights     = weights/sum(weights);

**/Users/florianoswald/actions-runner/_work/JPE-Ganapati-20241434/JPE-Ganapati-20241434/replication-package/Code/MonteCarlo/Script_QQ_simulation.m**

- Line 17, unix : addpath(project_path + '/functions/GFT_v2', '-frozen')
- Line 44, windows : fprintf('Running on %d cores, %d total economies\n', Ncores, Nsim);
- Line 141, unix : % Track converged/non-converged (only for the indices we ran)
- Line 148, windows : error('All iterations failed. First error:\n%s', err{find(~ok,1)});
- Line 193, windows : % fprintf(fid, 'i = %d\n%s\n\n', bad(t), err{bad(t)});
- Line 206, windows : fprintf('Completed in %.2f seconds. Successful: %d/%d\n', toc, numel(ok_idx), n_to_run);
- Line 268, windows : err_msg = sprintf('i=%d\n%s', i, getReport(ME,'extended','hyperlinks','off'));

**/Users/florianoswald/actions-runner/_work/JPE-Ganapati-20241434/JPE-Ganapati-20241434/replication-package/Code/MonteCarlo/solve_FP.m**

- Line 54, unix : N_i_sim = zeros(c, J); % firm entry/mass of new firms
- Line 121, windows : %fprintf('N_i: %.3e\n', any(isnan(N)));
- Line 140, windows : %fprintf('Error: %.3e\n', error_p);
- Line 159, windows : fprintf('Error: %.3e\n', error_w);
- Line 171, unix : % identify any non-finite along the key vectors/matrices
- Line 177, windows : fprintf('  offenders: X=%d N=%d P=%d w=%d\n', any(badX), any(badN), any(badP), any(badw));
- Line 206, windows : fprintf('[run %d] ERROR: %s - %s\n', b, ME.identifier, ME.message);
- Line 249, unix : % ---- run summary - identify converged/non-converged LN

**/Users/florianoswald/actions-runner/_work/JPE-Ganapati-20241434/JPE-Ganapati-20241434/replication-package/Code/p01_TRAINS_clean_data_v1_AAG.do**

- Line 38, unix : forvalues year=1995/2018 {
- Line 74, unix : forvalues year=1995/2018 {
- Line 167, unix : forvalues year=1995/2018 {
- Line 215, unix : forvalues year=1995/2018 {
- Line 266, unix : forvalues year=1995/2018 {
- Line 301, unix : forvalues year=1995/2018 {
- Line 340, unix : forvalues year=1995/2018 {

**/Users/florianoswald/actions-runner/_work/JPE-Ganapati-20241434/JPE-Ganapati-20241434/replication-package/Code/GMM_estimation_v3/GMM_gravity.m**

- Line 51, unix : kappa_r_implied = xi_tilde_epsilon_mean/xi_tilde_rho_mean;

**/Users/florianoswald/actions-runner/_work/JPE-Ganapati-20241434/JPE-Ganapati-20241434/replication-package/Code/p02_quantiles_colombia_firmexports.do**

- Line 13, unix : foreach p of numlist 1/99 {
- Line 20, unix : gen exp_sh = exp_agg/tot

**/Users/florianoswald/actions-runner/_work/JPE-Ganapati-20241434/JPE-Ganapati-20241434/replication-package/Code/ColombiaTests/Quantiles_Colombia_Wrapper.m**

- Line 52, unix : output_filename = '../../Output/QQ_Colombia_R2_Summary.csv';
- Line 56, windows : fprintf(fileID, '"%s",%f,%f,%f,%f\n', string(county_names{i, 1}), r2(i, 1), r2(i, 2), r2(i, 4), r2(i, 3));
- Line 58, windows : fprintf(fileID, '"Mean",%f,%f,%f,%f\n', mean_r2(1), mean_r2(2), mean_r2(4), mean_r2(3));

**/Users/florianoswald/actions-runner/_work/JPE-Ganapati-20241434/JPE-Ganapati-20241434/replication-package/Code/MonteCarlo/functions/inv_unique_interp.m**

- Line 18, windows : %fprintf('Computing: %d\n', k);

**/Users/florianoswald/actions-runner/_work/JPE-Ganapati-20241434/JPE-Ganapati-20241434/replication-package/Code/GMM_estimation_v3/p_overlay_elasticity.m**

- Line 66, unix : grad = (1/epsilon(j).^2) ;
- Line 75, unix : grad = [(rho(j)/epsilon(j).^2) ; (-1/epsilon(j))];

**/Users/florianoswald/actions-runner/_work/JPE-Ganapati-20241434/JPE-Ganapati-20241434/replication-package/Code/MonteCarlo/figures_gen.m**

- Line 481, windows : % fprintf('Spec %s: rmse_qq = %s, rmse_bl = %s\n', spec, mat2str(rmse_qq), mat2str(rmse_bl));
- Line 482, windows : % fprintf('  gft_qq_mse = %s, gft_bl_mse = %s\n', mat2str(spec_data.gft_qq_mse), mat2str(spec_data.gft_bl_mse));
- Line 483, windows : % fprintf('  mean(gft_qq_true) = %f, mean(gft_bl_true) = %f\n', ...
- Line 486, windows : fprintf('GFT Mean:\t%s: rmse_qq = %4.4f, rmse_bl = %4.4f\n', spec, mean(rmse_qq), mean(rmse_bl));
- Line 487, windows : fprintf('GFT Median:\t%s: rmse_qq = %4.4f, rmse_bl = %4.4f\n', spec, median(rmse_qq), median(rmse_bl));
- Line 562, unix : sgt = sgtitle(f, 'Gains from Trade: RMSE/Avg', 'FontWeight', 'bold', 'Interpreter', 'latex');

**/Users/florianoswald/actions-runner/_work/JPE-Ganapati-20241434/JPE-Ganapati-20241434/replication-package/Code/LogCorrected/plots.m**

- Line 8, unix : saveas(gcf,'../../Output/F1_plot_theta_e_vs_log_n','epsc')
- Line 11, unix : saveas(gcf,'../../Output/F1_plot_theta_c_vs_log_n','epsc')
- Line 14, unix : saveas(gcf,'../../Output/F1_plot_theta_vs_log_n','epsc')
- Line 17, unix : saveas(gcf,'../../Output/F1_plot_theta_i_vs_log_n','epsc')
- Line 21, unix : saveas(gcf,'../../Output/AF1_pdf_e','epsc')
- Line 24, unix : saveas(gcf,'../../Output/AF1_fEr','epsc')

**/Users/florianoswald/actions-runner/_work/JPE-Ganapati-20241434/JPE-Ganapati-20241434/replication-package/Code/Appendix_LogPareto/logpareto_inv_shifted.m**

- Line 15, unix : x = exp(k ./ (1 - p).^(1/alpha));

**/Users/florianoswald/actions-runner/_work/JPE-Ganapati-20241434/JPE-Ganapati-20241434/replication-package/Code/p04_Sample_Creation.do**

- Line 155, unix : replace BACI_v = BACI_v/1000
- Line 156, unix : replace BACI_q = BACI_q/1000
- Line 453, unix : replace CHN_X_ij_USD = CHN_X_ij_USD/1000000
- Line 466, unix : replace CHN_X_ij_USD = CHN_X_ij_yuan/US_FX if CHN_X_ij_USD == .
- Line 493, unix : replace AUS_X_ij = AUS_X_ij/US_FX
- Line 514, unix : global CEPII_Grav $ROOT/Data/CEPII/Gravity/gravdata_cepii.dta
- Line 520, unix : global CEPII_Grav $ROOT/Data/CEPII/Gravity_dta_V202102/Gravity_V202102.dta
- Line 742, unix : replace SurvivalRate`i'_combined = SurvivalRate`i'_combined/100
- Line 761, unix : gen N_ii_EDD_impute = EDD_N_iE/(export_probability/100)
- Line 800, unix : gen x_bar = TEC_X_ij/TEC_N_ij
- Line 808, unix : replace x_bar = A6i/1000000 if missing(x_bar)
- Line 811, unix : replace x_bar = AUS_X_ij/AUS_N_ij if missing(x_bar)
- Line 814, unix : replace x_bar = X_ij/AUS_N_ij if missing(x_bar)
- Line 821, unix : replace x_bar = X_ij/N_ij if iso3_o == iso3_d  & missing(x_bar)
- Line 830, unix : replace x_bar = X_ij/N_ij if missing(x_bar)
- Line 836, unix : gen n_ij_allsurvival = N_ij/N_ii
- Line 1001, unix : gen tariff_use = simpleAHS_w/100
- Line 1086, unix : gen tariff_use = simpleAHS_uw/100
- Line 1116, unix : gen tariff_use = TsimpleAHS_uw/100
- Line 1155, unix : gen tariff_use = TsimpleAHS_uw/100
- Line 1184, unix : gen tariff_use = simpleAHS_uw/100
- Line 1208, unix : gen ltariff_use = log(1+simpleAHS_uw/100)
- Line 1232, unix : gen ltariff_use = log(1+simpleAHS_uw/100)
- Line 1246, unix : gen ltariff_use = log(1+simpleAHS_uw/100)
- Line 1265, unix : gen tariff_use = simpleAHS_w/100
- Line 1280, unix : cap mkdir $ROOT/Data/Int/WIOD_sampleB/
- Line 1301, unix : cap mkdir $ROOT/Data/Int/WIOD_sampleB/
- Line 1302, unix : outsheet using $ROOT/Data/Int/WIOD_sampleB/${year}_TuwIV.csv, replace non nol comma
- Line 1343, unix : gen tariff_use = simpleAHS_uw/100
- Line 1486, unix : gen lGDP_diff = log(gdpcap_o/gdpcap_d)
- Line 1506, unix : gen X_ij_share = X_ij/X_ij_total
- Line 1507, unix : gen count_share = count/count_total
- Line 1561, unix : gen x_ij = X_ij/E_j
- Line 1564, unix : gen y_ij = X_ij/Y_i
- Line 1718, unix : gen kappa_i = Y_i/E_i

**/Users/florianoswald/actions-runner/_work/JPE-Ganapati-20241434/JPE-Ganapati-20241434/replication-package/Code/MasterProgram.do**

- Line 14, unix : global DATA   $ROOT/Data/
- Line 15, unix : global INT    $ROOT/Data/Int/
- Line 16, unix : global OUT    $ROOT/Output/
- Line 88, unix : do $ROOT/Code/p01_TRAINS_clean_data_v1_AAG.do
- Line 91, unix : do $ROOT/Code/p02_quantiles_colombia_firmexports.do
- Line 94, unix : do $ROOT/Code/p02_eora_trade_matrix_clean.do
- Line 97, unix : do $ROOT/Code/p02_ChinaData_Morrow.do
- Line 100, unix : do $ROOT/Code/p03_Create_Tariffs.do
- Line 103, unix : do $ROOT/Code/p04_Sample_Creation.do
- Line 104, unix : do $ROOT/Code/p05_Sample_Creation_HS2.do
- Line 107, unix : do $ROOT/Code/p06_ReducedForm.do
- Line 112, unix : cd $ROOT/Code/GMM_estimation_v3/
- Line 116, unix : cd $ROOT/Code/GFT_v2/
- Line 120, unix : cd $ROOT/Code/LinearCF_v3/
- Line 125, unix : cd $ROOT/Code/ColombiaTests/
- Line 129, unix : cd $ROOT/Code/MonteCarlo/
- Line 133, unix : cd $ROOT/Code/LogCorrected
- Line 137, unix : cd $ROOT/Code/Appendix_LogPareto

**/Users/florianoswald/actions-runner/_work/JPE-Ganapati-20241434/JPE-Ganapati-20241434/replication-package/Code/GMM_estimation_v3/graph_overlay.m**

- Line 47, unix : grad = [(1/epsilon1(j).^2) ];
- Line 57, unix : grad = [(rho1(j)/epsilon1(j).^2) ; (-1/epsilon1(j))];
- Line 286, unix : case  'output/F2ab_dual_'
- Line 288, unix : case  'output/F4b_HS_origin'

**/Users/florianoswald/actions-runner/_work/JPE-Ganapati-20241434/JPE-Ganapati-20241434/replication-package/Code/p03_Create_Tariffs.do**

- Line 89, unix : forvalues year=2012/2012 {
- Line 124, unix : forvalues year=2012/2012 {
- Line 137, unix : gen hs2 = floor(hs92/1000)

**/Users/florianoswald/actions-runner/_work/JPE-Ganapati-20241434/JPE-Ganapati-20241434/replication-package/Code/p06_ReducedForm.do**

- Line 13, unix : replace rank_n_ij = rank_n_ij/_N
- Line 17, unix : gen tariff_use = TsimpleAHS_uw/100
- Line 42, unix : replace rank_n_ij = rank_n_ij/NN

**/Users/florianoswald/actions-runner/_work/JPE-Ganapati-20241434/JPE-Ganapati-20241434/replication-package/Code/GMM_estimation_v3/GenerateBaseline.m**

- Line 5, unix : addpath('../Estimation/GMM_estimation_v3')
- Line 6, unix : data_dir = '../../Data/Int/WIOD_sampleB/';

**/Users/florianoswald/actions-runner/_work/JPE-Ganapati-20241434/JPE-Ganapati-20241434/replication-package/Code/GFT_v2/GFT_Script_v17_EKK.m**

- Line 20, unix : z           = unif_EKK.^(-1/theta_EKK);
- Line 32, unix : elast_r_EKK =     (r_vec/iter)./cumtrapz(unif_EKK,(r_vec/iter))./unif_EKK-1;
- Line 33, unix : e               = e_vec/iter;
- Line 35, unix : r               = r_vec/iter;

**/Users/florianoswald/actions-runner/_work/JPE-Ganapati-20241434/JPE-Ganapati-20241434/replication-package/Code/GMM_estimation_v3/make_data_combo.m**

- Line 251, unix : disp('discrete measure of China (1/0)')

**/Users/florianoswald/actions-runner/_work/JPE-Ganapati-20241434/JPE-Ganapati-20241434/replication-package/Code/LinearCF_v3/LinearCF.m**

- Line 5, unix : data_dir = [ '../../Data/Int/WIOD_sampleB/'];
- Line 104, unix : load('../GMM_estimation_v3/gamma_gravity.mat');
- Line 138, unix : G_ijr        = csvread([  '../../Data/Int/WIOD_sampleB/G_ij_'     file_base ]);
- Line 152, unix : G_ijr        = csvread([  '../../Data/Int/WIOD_sampleB/G_ij_'     file_base ]);
- Line 169, unix : G_ijr        = csvread([  '../../Data/Int/WIOD_sampleB/G_ij_'     file_base ]);

**/Users/florianoswald/actions-runner/_work/JPE-Ganapati-20241434/JPE-Ganapati-20241434/replication-package/Code/LogCorrected/plot_logcorr_pareto.m**

- Line 78, unix : tU   = max(tL + 1, 1 - ln_p/a + 10);

**/Users/florianoswald/actions-runner/_work/JPE-Ganapati-20241434/JPE-Ganapati-20241434/replication-package/Code/LogCorrected/plot_theta_i_vs_log_n.m**

- Line 17, unix : %   theta^e(n) = alpha^e + gamma^e/ln(epsilon(n))

**/Users/florianoswald/actions-runner/_work/JPE-Ganapati-20241434/JPE-Ganapati-20241434/replication-package/Code/MonteCarlo/Script_GMM_simulation.m**

- Line 17, unix : addpath(project_path + '/functions/GFT_v2', '-frozen')
- Line 44, windows : fprintf('Running on %d cores, %d total economies\n', Ncores, Nsim);
- Line 139, unix : % Track converged/non-converged (only for the indices we ran)
- Line 146, windows : error('All iterations failed. First error:\n%s', err{find(~ok,1)});
- Line 190, windows : printf(fid, 'i = %d\n%s\n\n', bad(t), err{bad(t)});
- Line 201, windows : fprintf('Completed in %.2f seconds. Successful: %d/%d\n', toc, numel(ok_idx), n_to_run);
- Line 269, windows : err_msg = sprintf('i=%d\n%s', i, getReport(ME,'extended','hyperlinks','off'));

**/Users/florianoswald/actions-runner/_work/JPE-Ganapati-20241434/JPE-Ganapati-20241434/replication-package/Code/GFT_v2/GFT_script_v17.m**

- Line 7, unix : data_dir  =  '../../Data/Int/WIOD_sampleB/';
- Line 8, unix : excel_output = [ '../../Output/GFT_results_revision_' year '.xlsx'];
- Line 28, unix : load( '../GMM_estimation_v3/gamma_gravity.mat');
- Line 86, windows : fprintf( 'ERROR: % 10.5f\t\n',    mean(abs(W_hat_gen_pareto./W_hat_pareto_simple))-1)
- Line 87, windows : fprintf( 'ERROR: % 10.5f\t\n',    mean(abs(W_hat_pareto_simple./W_hat_gen_pareto_alt))-1)
- Line 153, unix : weights     = weights/sum(weights);
- Line 191, unix : xlabel('Import Share $\log(1/x_{ii})$','FontSize',16,'Interpreter','Latex');
- Line 332, unix : xlabel('Import Share $\log(1/x_{ii})$','FontSize',16,'Interpreter','Latex');

**/Users/florianoswald/actions-runner/_work/JPE-Ganapati-20241434/JPE-Ganapati-20241434/replication-package/Code/MonteCarlo/generate_QQ_plot.m**

- Line 26, windows : fprintf('Simulations to plot: %s\n', mat2str(sim_indices));
- Line 31, windows : fprintf('\n\n');
- Line 60, windows : fprintf('Plot naming convention: QQ_scatter_[fct_form]_[basis]_k[knots]_sim[id]_i[ii]_j[jj].pdf\n');
- Line 68, unix : % estimation_path = 'path/to/estimation';
- Line 69, unix : % data_path = 'path/to/data';
- Line 70, unix : % output_path = 'path/to/plots';
- Line 98, unix : % load('path/to/Simulated_data_LogNormal_alt.mat', 'M', 'Q_lnx_all');
- Line 111, unix : % % Check output: output/QQ_scatter_LogNormal_LogNormal_k1_sim1_i1_j2.pdf

**/Users/florianoswald/actions-runner/_work/JPE-Ganapati-20241434/JPE-Ganapati-20241434/replication-package/Code/Appendix_LogPareto/AppendixFigure_LogPareto.m**

- Line 43, unix : saveas(gca,'../../Output/FA_logPareto_Sim','epsc')

