## Potentially Hardcoded Numeric Constants


We found the following set of hard coded numbers. This may be completely legitimate (parameter input, thresholds for computations, etc), and is hence only for information.

**/Users/florianoswald/actions-runner/_work/JPE-Ganapati-20241434/JPE-Ganapati-20241434/replication-package/Code/MonteCarlo/get_fct_form_params.m**

- Line 27, : params.sigma_r = 0.000;
- Line 28, : params.sigma_f = 0.000;

**/Users/florianoswald/actions-runner/_work/JPE-Ganapati-20241434/JPE-Ganapati-20241434/replication-package/Code/GMM_estimation_v3/Master_script_gravity_v3.m**

- Line 40, : %    (from Teti Tariffs) gamma_gravity = [ 0.3520   -0.0214   -0.0189   -0.2324   -0.1502    0.2812   -0.0834   -0.0015   -0.0576   -0.0471];
- Line 41, : %    (from Base Tariffs) gamma_gravity = [ 0.3920   -0.0408   -0.0131   -0.2596   -0.1635    0.3402   -0.1012    0.0039   -0.0985   -0.0671];
- Line 42, : %    (from Replication ) gamma_gravity = [ 0.3558   -0.0216   -0.0200   -0.2365   -0.1511    0.2860   -0.0832   -0.0047   -0.0636   -0.0478];
- Line 243, : % gamma_alt_UW =    0.3161    0.2699
- Line 490, : gamma_guessHS = [0.5533   -0.0410   -0.0557   -0.4006   -0.3545    0.2085   -0.0830   -0.0148   -0.0016    0.06080];

**/Users/florianoswald/actions-runner/_work/JPE-Ganapati-20241434/JPE-Ganapati-20241434/replication-package/Code/GMM_estimation_v3/CC_se_thetae.m**

- Line 7, : percentile_se = norminv([0.025 0.975]);
- Line 13, : alpha = [ 0.9500 ];

**/Users/florianoswald/actions-runner/_work/JPE-Ganapati-20241434/JPE-Ganapati-20241434/replication-package/Code/MonteCarlo/functions/figure_plotter_QQ_minimal.m**

- Line 36, : set(ax, 'Units', 'normalized', 'Position', [0.13 0.11 0.775 0.815], ...

**/Users/florianoswald/actions-runner/_work/JPE-Ganapati-20241434/JPE-Ganapati-20241434/replication-package/Code/GFT_v2/GFT_division2.m**

- Line 3, : a  = scatter(xmetric(rich),ymetric1(rich),[],[0.6350 0.0780 0.1840],'filled','o');
- Line 6, : b  = scatter(xmetric(poor),ymetric1(poor),[],[0 0.4470 0.7410],'filled','square');

**/Users/florianoswald/actions-runner/_work/JPE-Ganapati-20241434/JPE-Ganapati-20241434/replication-package/Code/LogCorrected/plot_theta_e_vs_log_n.m**

- Line 15, : n_min = 0.0001;

**/Users/florianoswald/actions-runner/_work/JPE-Ganapati-20241434/JPE-Ganapati-20241434/replication-package/Code/GMM_estimation_v3/CC_se_splines.m**

- Line 9, : xtemp = (0.05:0.001:0.95)';
- Line 10, : %         xtemp = (0.01:0.001:0.99)';
- Line 15, : %         scale_upper = quantile(d.R_nE_ij, 0.999);
- Line 16, : %         scale_lower = quantile(d.R_nE_ij, 0.001);
- Line 35, : alpha = [ 0.9500 ];
- Line 88, : percentile_se = norminv([0.025 0.975]);

**/Users/florianoswald/actions-runner/_work/JPE-Ganapati-20241434/JPE-Ganapati-20241434/replication-package/Code/GMM_estimation_v3/make_elasticity_combo.m**

- Line 75, : percentile_se = norminv([0.025 0.975]);

**/Users/florianoswald/actions-runner/_work/JPE-Ganapati-20241434/JPE-Ganapati-20241434/replication-package/Code/LogCorrected/plot_theta_combined.m**

- Line 21, : n_min = 0.0001;

**/Users/florianoswald/actions-runner/_work/JPE-Ganapati-20241434/JPE-Ganapati-20241434/replication-package/Code/GMM_estimation_v3/Simulation_Data_v10.m**

- Line 11, : tau             = 1.83; % https://pubs.aeaweb.org/doi/pdfplus/10.1257/aer.20130351 page 1128
- Line 59, : % That paper doesn't estimate it, it simply uses Head Meyer AER where 0.797
- Line 109, : H                  = 2.85; % From https://pubs.aeaweb.org/doi/pdfplus/10.1257/aer.20130351
- Line 178, : ax.XTick   = [-9.210340372 -6.907755279 -5.298317367 -4.605170186 -2.995732274 -2.302585093 -0.693147181 0 ];
- Line 196, : ax.XTick   = [-9.210340372 -6.907755279 -5.298317367 -4.605170186 -2.995732274 -2.302585093 -0.693147181 0 ];
- Line 216, : ax.XTick   = [-9.210340372 -6.907755279 -5.298317367 -4.605170186 -2.995732274 -2.302585093 -0.693147181 0 ];
- Line 236, : ax.XTick   = [-9.210340372 -6.907755279 -5.298317367 -4.605170186 -2.995732274 -2.302585093 -0.693147181 0 ];
- Line 255, : ax.XTick   = [-9.210340372 -6.907755279 -5.298317367 -4.605170186 -2.995732274 -2.302585093 -0.693147181 0 ];
- Line 274, : ax.XTick   = [-9.210340372 -6.907755279 -5.298317367 -4.605170186 -2.995732274 -2.302585093 -0.693147181 0 ];
- Line 348, : ax.XTick   = [-9.210340372 -6.907755279 -5.298317367 -4.605170186 -2.995732274 -2.302585093 -0.693147181 0 ];
- Line 369, : ax.XTick   = [-9.210340372 -6.907755279 -5.298317367 -4.605170186 -2.995732274 -2.302585093 -0.693147181 0 ];
- Line 389, : ax.XTick   = [-9.210340372 -6.907755279 -5.298317367 -4.605170186 -2.995732274 -2.302585093 -0.693147181 0 ];
- Line 411, : ax.XTick   = [-9.210340372 -6.907755279 -5.298317367 -4.605170186 -2.995732274 -2.302585093 -0.693147181 0 ];
- Line 432, : ax.XTick   = [-9.210340372 -6.907755279 -5.298317367 -4.605170186 -2.995732274 -2.302585093 -0.693147181 0 ];
- Line 459, : ax.XTick   = [-9.210340372 -6.907755279 -5.298317367 -4.605170186 -2.995732274 -2.302585093 -0.693147181 0 ];

**/Users/florianoswald/actions-runner/_work/JPE-Ganapati-20241434/JPE-Ganapati-20241434/replication-package/Code/MonteCarlo/solve_FP.m**

- Line 42, : n_ij        = @(fbar, rbar, w_i, w_j, P_j) max(min(1 - H_e(log((sigma .* fbar .* w_i.^sigma)./(rbar .* P_j .* w_j))), 0.999), 1e-18);

**/Users/florianoswald/actions-runner/_work/JPE-Ganapati-20241434/JPE-Ganapati-20241434/replication-package/Code/GFT_v2/GFT_division.m**

- Line 3, : a  = scatter(xmetric(rich),ymetric1(rich),[],[0.6350 0.0780 0.1840],'filled','o');
- Line 6, : b  = scatter(xmetric(poor),ymetric1(poor),[],[0 0.4470 0.7410],'filled','square');

**/Users/florianoswald/actions-runner/_work/JPE-Ganapati-20241434/JPE-Ganapati-20241434/replication-package/Code/GMM_estimation_v3/CC_se.m**

- Line 17, : alpha = [ 0.9500 ];
- Line 71, : percentile_se = norminv([0.025 0.975]);

**/Users/florianoswald/actions-runner/_work/JPE-Ganapati-20241434/JPE-Ganapati-20241434/replication-package/Code/MonteCarlo/functions/GFT_fct_general_fast.m**

- Line 33, : n_grid = logspace(log(min(min(n_ij))), 0, 500); % Log-spaced from 0.0001 to 1

**/Users/florianoswald/actions-runner/_work/JPE-Ganapati-20241434/JPE-Ganapati-20241434/replication-package/Code/MonteCarlo/functions/evknots.m**

- Line 59, : % k=[-9.88055994 -8	-7.210481 -5	-3.66798696 ]; % hardcoded from the paper
- Line 60, : % k=[-9.88055994 -8.5	-7.810481 -5	-3.66798696 ]; % hardcoded from the paper
- Line 61, : % [-9.8806 -7.0210 -3.6680]
- Line 84, : k=[-9.88055994	-7.0210481	-3.66798696 ]; % hardcoded from the paper

**/Users/florianoswald/actions-runner/_work/JPE-Ganapati-20241434/JPE-Ganapati-20241434/replication-package/Code/GMM_estimation_v3/p_overlay_elasticity.m**

- Line 38, : % percentile_se = norminv([0.025 0.975]);
- Line 448, : ax.XTick   = [-9.210340372 -6.907755279  -4.605170186  -2.302585093  0 ];

**/Users/florianoswald/actions-runner/_work/JPE-Ganapati-20241434/JPE-Ganapati-20241434/replication-package/Code/MonteCarlo/figures_gen.m**

- Line 254, : xt = log(0.0001):log(10):0;
- Line 358, : ax.XTick   = [-9.210340372 -6.907755279  -4.605170186  -2.302585093  0 ];
- Line 689, : set(newAx, 'Units', 'normalized', 'Position', [0.13 0.11 0.775 0.815], ...

**/Users/florianoswald/actions-runner/_work/JPE-Ganapati-20241434/JPE-Ganapati-20241434/replication-package/Code/LinearCF_v3/CF_agg_plot.m**

- Line 3, : a = scatter(ln_n_ij(group1),metric(group1),[],[0.6350 0.0780 0.1840],'filled','o');
- Line 5, : b = scatter(ln_n_ij(group2),metric(group2),[],[0 0.4470 0.7410],'filled','square');

**/Users/florianoswald/actions-runner/_work/JPE-Ganapati-20241434/JPE-Ganapati-20241434/replication-package/Code/p04_Sample_Creation.do**

- Line 1670, : replace G4_ij = 3 if rich_orig == 0 & gdpcap_o <= 2995 // 7.3140e+03
- Line 1671, : replace G4_ij = 2 if rich_orig == 0 & gdpcap_o <= 755 // 2.9217e+03

**/Users/florianoswald/actions-runner/_work/JPE-Ganapati-20241434/JPE-Ganapati-20241434/replication-package/Code/GMM_estimation_v3/graph_overlay.m**

- Line 21, : %         percentile_se = norminv([0.025 0.975]);
- Line 120, : %                 a(num+1)=plot(ln_ij,elast_e_pareto(range1),'-','LineWidth',2,'Color',[0 0.4470 0.7410]	);
- Line 121, : a(num+1)=plot(ln_ij,elast_e_tpareto(range1),'--','LineWidth',2,'Color',[0.8500 0.3250 0.0980]	);
- Line 122, : a(num+2)=plot(ln_ij,elast_e_logn(range1),':','LineWidth',2,'Color',[0.9290 0.6940 0.1250]);
- Line 123, : a(num+3)=plot(log(unif_EKK(range2)),smoothdata(elast_e_EKK(range2),'gaussian',20),'-.','LineWidth',2,'Color',[0.4940 0.1840 0.5560]);
- Line 139, : %                 a(num+1)=plot(ln_ij,elast_r_pareto(range1),'-','LineWidth',2,'Color',[0 0.4470 0.7410]	);
- Line 140, : a(num+1)=plot(ln_ij,elast_r_tpareto(range1),'--','LineWidth',2,'Color',[0.8500 0.3250 0.0980]	);
- Line 141, : a(num+2)=plot(ln_ij,elast_r_logn(range1),':','LineWidth',2,'Color',[0.9290 0.6940 0.1250]);
- Line 142, : a(num+3)=plot(log(unif_EKK(range2)),smoothdata(elast_r_EKK(range2),'gaussian',20),'-.','LineWidth',2,'Color',[0.4940 0.1840 0.5560]	);
- Line 161, : %                 a(num+1)=plot(ln_ij,theta_pareto(range1),'-','LineWidth',2,'Color',[0 0.4470 0.7410]	);
- Line 162, : a(num+1)=plot(ln_ij,theta_tpareto(range1),'--','LineWidth',2,'Color',[0.8500 0.3250 0.0980]	);
- Line 163, : a(num+2)=plot(ln_ij,theta_logn(range1),':','LineWidth',2,'Color',[0.9290 0.6940 0.1250]);
- Line 164, : a(num+3)=plot(log(unif_EKK(range2)),smoothdata(theta_EKK(range2),'gaussian',20),'-.','LineWidth',2,'Color',[0.4940 0.1840 0.5560]		);
- Line 251, : %             a(num+1)=plot(ln_ij,(1-d(i).sigma)*1./elast_e_pareto(range1),'-','LineWidth',2,'Color',[0 0.4470 0.7410]);
- Line 252, : a(num+1)=plot(ln_ij,(1-d(i).sigma)*1./elast_e_tpareto(range1),'-','LineWidth',2,'Color',[0.8500 0.3250 0.0980]);
- Line 253, : a(num+2)=plot(ln_ij,(1-d(i).sigma)*1./elast_e_logn(range1),'-','LineWidth',2,'Color',[0.9290 0.6940 0.1250]);
- Line 254, : a(num+3)=plot(log(unif_EKK(range2)),(1-d(i).sigma)*1./smoothdata(elast_e_EKK(range2),'gaussian',10),'-','LineWidth',2,'Color',[0.4940 0.1840 0.5560]	);
- Line 269, : ax.XTick   = [-9.210340372 -6.907755279  -4.605170186  -2.302585093  0 ];

**/Users/florianoswald/actions-runner/_work/JPE-Ganapati-20241434/JPE-Ganapati-20241434/replication-package/Code/LogCorrected/plot_theta_c_vs_log_n.m**

- Line 14, : n_min = 0.0001;

**/Users/florianoswald/actions-runner/_work/JPE-Ganapati-20241434/JPE-Ganapati-20241434/replication-package/Code/LinearCF_v3/CF_diagnostic.m**

- Line 2, : a = scatter(ln_n_ij(rich),plot_var(rich),[],[0.6350 0.0780 0.1840],'filled','o');
- Line 4, : b = scatter(ln_n_ij(poor),plot_var(poor),[],[0 0.4470 0.7410],'filled','square');

**/Users/florianoswald/actions-runner/_work/JPE-Ganapati-20241434/JPE-Ganapati-20241434/replication-package/Code/p06_ReducedForm.do**

- Line 21, : gen dist_med = cond(ldistw> 8.562,1,0)

**/Users/florianoswald/actions-runner/_work/JPE-Ganapati-20241434/JPE-Ganapati-20241434/replication-package/Code/GFT_v2/GFT_agg_plot.m**

- Line 3, : a = scatter(ln_n_ij(rich),metric(rich),[],[0.6350 0.0780 0.1840],'filled','o');
- Line 5, : b = scatter(ln_n_ij(poor),metric(poor),[],[0 0.4470 0.7410],'filled','square');

**/Users/florianoswald/actions-runner/_work/JPE-Ganapati-20241434/JPE-Ganapati-20241434/replication-package/Code/GFT_v2/GFT_Script_v17_EKK.m**

- Line 83, : a = scatter(ln_n_ij(rich),metric(rich),[],[0.6350 0.0780 0.1840],'filled','o');
- Line 85, : b = scatter(ln_n_ij(poor),metric(poor),[],[0 0.4470 0.7410],'filled','square');

**/Users/florianoswald/actions-runner/_work/JPE-Ganapati-20241434/JPE-Ganapati-20241434/replication-package/Code/GMM_estimation_v3/GenerateBaseline.m**

- Line 11, : % gamma_gravity =    [0.3520   -0.0214   -0.0189   -0.2324   -0.1502    0.2812   -0.0834   -0.0015   -0.0576   -0.0471];

**/Users/florianoswald/actions-runner/_work/JPE-Ganapati-20241434/JPE-Ganapati-20241434/replication-package/Code/LinearCF_v3/LinearCF.m**

- Line 107, : %    (from Teti Tariffs) gamma_gravity = [ 0.3520   -0.0214   -0.0189   -0.2324   -0.1502    0.2812   -0.0834   -0.0015   -0.0576   -0.0471];
- Line 108, : %    (from Base Tariffs) gamma_gravity = [ 0.3920   -0.0408   -0.0131   -0.2596   -0.1635    0.3402   -0.1012    0.0039   -0.0985   -0.0671];

**/Users/florianoswald/actions-runner/_work/JPE-Ganapati-20241434/JPE-Ganapati-20241434/replication-package/Code/LogCorrected/plot_logcorr_pareto.m**

- Line 33, : n  = linspace(1e-4, 0.9999, 1500);

**/Users/florianoswald/actions-runner/_work/JPE-Ganapati-20241434/JPE-Ganapati-20241434/replication-package/Code/LogCorrected/plot_theta_i_vs_log_n.m**

- Line 21, : n_min = 0.0001;

**/Users/florianoswald/actions-runner/_work/JPE-Ganapati-20241434/JPE-Ganapati-20241434/replication-package/Code/GFT_v2/GFT_script_v17.m**

- Line 26, : % gamma_gravity = [ 0.3920   -0.0408   -0.0131   -0.2596   -0.1635    0.3402   -0.1012    0.0039   -0.0985   -0.0671];
- Line 27, : % gamma_gravity =    [0.3520   -0.0214   -0.0189   -0.2324   -0.1502    0.2812   -0.0834   -0.0015   -0.0576   -0.0471];
- Line 436, : tau        = 1.83; % https://pubs.aeaweb.org/doi/pdfplus/10.1257/aer.20130351 page 1128
- Line 440, : H          = 2.85; % From https://pubs.aeaweb.org/doi/pdfplus/10.1257/aer.20130351
- Line 502, : % https://pubs.aeaweb.org/doi/pdfplus/10.1257/aer.104.5.310
- Line 593, : a  = scatter(xmetric(rich),ymetric1(rich),[],[0.6350 0.0780 0.1840],'filled','o');
- Line 596, : b  = scatter(xmetric(poor),ymetric1(poor),[],[0 0.4470 0.7410],'filled','square');
- Line 607, : a  = scatter(xmetric(rich),ymetric1(rich),[],[0.6350 0.0780 0.1840],'filled','o');
- Line 610, : b  = scatter(xmetric(poor),ymetric1(poor),[],[0 0.4470 0.7410],'filled','square');
- Line 621, : a  = scatter(xmetric(rich),ymetric1(rich),[],[0.6350 0.0780 0.1840],'filled','o');
- Line 624, : b  = scatter(xmetric(poor),ymetric1(poor),[],[0 0.4470 0.7410],'filled','square');
- Line 639, : a  = scatter(xmetric(rich),ymetric1(rich),[],[0.6350 0.0780 0.1840],'filled','o');
- Line 642, : b  = scatter(xmetric(poor),ymetric1(poor),[],[0 0.4470 0.7410],'filled','square');
- Line 659, : a  = scatter(xmetric(rich),ymetric1(rich),[],[0.6350 0.0780 0.1840],'filled','o');
- Line 662, : b  = scatter(xmetric(poor),ymetric1(poor),[],[0 0.4470 0.7410],'filled','square');

**/Users/florianoswald/actions-runner/_work/JPE-Ganapati-20241434/JPE-Ganapati-20241434/replication-package/Code/LogCorrected/plot_Er_vs_log_n.m**

- Line 13, : n_min = 0.0001;

**/Users/florianoswald/actions-runner/_work/JPE-Ganapati-20241434/JPE-Ganapati-20241434/replication-package/Code/LinearCF_v3/CF_flow_plot.m**

- Line 4, : scatter(n(IND==0),(plot_var(IND==0)),'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2,'MarkerEdgeColor',C1,'MarkerFaceColor',C1);   a = scatter([],[],'MarkerEdgeColor',C1,'MarkerFaceColor',C1);%a.Color = C1;%[0 0.4470 0.7410];
- Line 5, : scatter(n(IND==1),(plot_var(IND==1)),'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2,'MarkerEdgeColor',C1);                        b = scatter([],[],'MarkerEdgeColor',C2);%b.Color = C2;%[0.8500 0.3250 0.0980];
- Line 6, : scatter(n(IND==2),(plot_var(IND==2)),'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2,'MarkerEdgeColor',C3,'MarkerFaceColor',C3);   c = scatter([],[],'MarkerEdgeColor',C3,'MarkerFaceColor',C3);%c.Color = C3;%[0.9290 0.6940 0.1250];
- Line 7, : scatter(n(IND==3),(plot_var(IND==3)),'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2,'MarkerEdgeColor',C4);                        d = scatter([],[],'MarkerEdgeColor',C4);%d.Color = C4;%[0.4940 0.1840 0.5560];
- Line 11, : f =gca; f.XTick   = [-9.210340372 -6.907755279  -4.605170186  -2.302585093  0 ];f.XTickLabel = {'0.01%' ,'0.1%','1%','10%','100%'};

**/Users/florianoswald/actions-runner/_work/JPE-Ganapati-20241434/JPE-Ganapati-20241434/replication-package/Code/Appendix_LogPareto/AppendixFigure_LogPareto.m**

- Line 41, : ax.XTick   = [-9.210340372 -6.907755279  -4.605170186  -2.302585093  0 ];

