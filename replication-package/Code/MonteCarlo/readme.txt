Monte Carlo Simulation Code 

This readme describes the master script and the main scripts it calls for running Monte Carlo simulations to estimate trade elasticities using GMM and QQ methods across different functional form specifications. I used Claude Opus to help generate this readme file.

================================================================================
MASTER FILE
================================================================================

1. run_mc_master.m

Description:
Master script that executes the workflow: data generation, 
GMM estimation, QQ estimation, plotting. Supports four functional form specifications: 
'LogNormal' (decreasing elasticity), 
'MPareto' (modified Pareto), and 'MEstimates' (increasing elasticity).

Input:
- User-specific project path and computational settings (cores, number of simulations)
- Specification-estimator mapping structures defining which estimators to run
- flags: create_simulation, run_gmm, run_qq, do_GFT, save_graph

Output:
- Simulated data files for each specification
- GMM estimation results (trade elasticities, extensive and intensive margin elasticities/functions, GFT)
- QQ estimation results (trade elasticities, extensive and intensive margin elasticities/functions, GFT)
- Plots (elasticity estimates, GFT histograms, QQ plots)
- All outputs saved in /output/data/, /output/estimation/, and /output/plots/

================================================================================
DATA GENERATION SCRIPTS (Called in the specified order)
================================================================================

1. get_fct_form_params.m

Description:
Returns default parameter values for each functional form specification.

Input:
- fct_form: String specifying functional form ('LogNormal', 'MPareto', 'MEstimates')

Output:
- params: Struct containing all model parameters (sigma, nu_ij, mu_f, sigma_z, etc.)


2. create_elasticity_functions.m

Description:
Generates elasticity functions (epsilon, rho, kappa) and other equilibrium objects as function 
handles based on the specified functional form and parameters.

Input:
- fct_form: Functional form specification
- project_path: Path to project directory
- output_path: Path to output directory
- params: Structure from get_fct_form_params.m

Output:
- out: Struct containing function handles for:
  - epsilon: extensive margin elasticity function
  - eps_true: true extensive margin elasticity values
  - rho: intensive margin elasticity function
  - r_true: true intensive margin elasticity
  - kappa: integral of the ratio of epsilon to rho 
  - H_e: CDF of entry potentials
  - Q_lnx: Sales distribution quantile function
  - p_grid: Percentile grid (1-99%)


3. solve_FP.m

Description:
Solves for general equilibrium via fixed point iteration. Generates 
simulated trade data

Input:
- project_path, output_path: Directory paths
- data_params: Parameter structure from get_fct_form_params.m
- epsilon, eps_true, rho, r_true, kappa, H_e, Q_lnx: Function handles from create_elasticity_functions.m
- p_grid: Percentile grid
- sims: Number of Monte Carlo simulations
- cores: Number of cores for parallelization
- upd_step (gamma): Updating parameter for fixed point iteration (typically 0.1)

Output:
- Simulated_data_[fct_form]_alt.mat containing:
  - M: Matrix with bilateral trade variables (c^2 × 7 × sims)
    Columns: [i_idx, j_idx, R_n_ij, R_xbar_ij, Z1, Z2, X_ij]
  - w_i, P_j, N_i (equilibrium wages, prices, mass of firms)
  - epsilon_true, rho_true: extensive/intensive margin elasticity functions (stored as handles)
  - epsilon_true_el, rho_true_el, theta_true: extensive/intensive margin and trade elasticities over equilibrium support
  - Q_lnx_all: Sales distribution quantiles (c^2 × 99 × sims)
  - conv_flag, nan_flag: Convergence diagnostics

================================================================================
ESTIMATION SCRIPTS
================================================================================

1. Script_GMM_simulation.m

Description:
Performs GMM estimation of trade elasticities. 
Supports different basis functions (Spline, LogNormal) and varying numbers of knots. 
Can optionally compute Gains from Trade (GFT). Supports running specific simulations 
or ranges via 'RunSingle' and 'RunRange' options.

Input:
- project_path, output_path: Directory paths
- fct_form: Functional form specification
- estimators: Cell array of estimator structures with fields {basis, knots}
- cores: Number of cores for parallel processing
- sims: Total number of simulations available
- do_GFT: Flag to compute gains from trade (0 or 1)
- Optional: 'RunSingle', sim_id or 'RunRange', [start, end]

Output:
- Simulated_[fct_form]_Estimates_[basis]_[knots]_alt.mat containing:
  - Estimated (epsilon_n, rho_n, theta_n) and true (epsilon_true_el, rho_true_el, theta_true) elasticities in struct p
  - Estimated parameters and functions for each simulation and elasticity functions in struct o
  - Rn_mat: Log entry probabilities used in estimation
  - Optional: estimated (w_hat) and true (w_hat_true) GFT if do_GFT=1


2. Script_QQ_simulation.m

Description:
Performs QQ (Quantile-Quantile) regression estimation using sales distribution quantiles. 
Always uses LogNormal basis with 1 knot to recover elasticities. The 'force_n' parameter determines whether to estimate with fixed n=1 assumption ('n1') or variable n ('nv'). Can optionally compute GFT.

Input:
- project_path, output_path: Directory paths
- fct_form: Functional form specification
- estimators: Cell array with fields {basis, knots, force_n}
  - basis: Must be 'LogNormal'
  - knots: Must be 1
  - force_n: 'n1' (assume n=1) or 'nv' (estimate n)
- cores: Number of cores for parallel processing
- sims: Total number of simulations available
- doGFT: Flag to compute gains from trade
- Optional: 'RunSingle', sim_id or 'RunRange', [start, end]

Output:
- Simulated_[fct_form]_Estimates_QQLogNormal_1_[force_n]_alt.mat containing:
  - Estimated (epsilon_n, rho_n, theta_n) and true (epsilon_true_el, rho_true_el, theta_true) elasticities in struct p
  - Estimated dispersion parameter in struct o
  - R^2 from minimal distance estimator in struct p
  - Rn_mat: Entry probabilities
  - Optional: estimated (w_hat) and true (w_hat_true) GFT if do_GFT=1

================================================================================
PLOTTING SCRIPTS
================================================================================

1. figures_gen.m

Description:
Creates elasticity plots comparing GMM and QQ estimates against true 
values, with support for multiple specifications and estimators. Optionally creates GFT histograms.

Input:
- baseline_path: Path to estimation results (/output/estimation)
- output_path: Path to save plots (/output/plots)
- version: Version string (must be "alt")
- qq: Structure mapping specifications to QQ force_n options (e.g., struct('LogNormal','nv'))
- Optional name-value pairs:
  - 'specifications': Cell array of specs to plot (default: {'LogNormal', 'MEstimates'})
  - 'gmm_estimators': Cell array of GMM estimator structures to include
  - 'save_separate': Save individual plots (default: true)
  - 'save_combined': Save combined comparison plot (default: true)
  - 'plot_gft': Include GFT histograms (default: false)
  - 'ylim_theta': Y-axis limits for theta plot (default: [0 8])
  - 'xlim_gft': X-axis limits for GFT plot (default: [0 8])
  - 'font_size': Base font size (default: 12)

Output:
- Individual elasticity plots: Elasticity_[spec]_[suffix].pdf
- GFT histograms: GFT_[spec].pdf (if plot_gft=true)


2. generate_QQ_plot.m

Description:
Generates detailed QQ scatter plots for specified simulations and country pairs, plotting simulated and estimated quantiles of the sales distribution.

Input:
- estimation_path: Path to estimation results (can be empty, not currently used)
- data_path: Path to simulated data (/output/data)
- output_path: Path to save plots
- fct_form: Functional form specification
- basis: Basis function (must be 'LogNormal' for QQ)
- knots: Number of knots (typically 1)
- force_n: QQ option ('n1' or 'nv')
- sim_indices: Vector of simulation IDs to plot (e.g., [1, 5])
- plot_pairs: Cell array of [i,j] country pairs to plot (e.g., {[1,4]})

Output:
- QQ scatter plots: QQ_scatter_[fct_form]_[basis]_k[knots]_sim[id]_i[ii]_j[jj].pdf
  Each plot shows simulated quantiles vs estimated quantiles for the specified 
  country pair and simulation.

================================================================================
WORKFLOW SUMMARY
================================================================================

The typical workflow executed by run_mc_master.m:

1. DATA GENERATION (if create_simulation=1):
   For each functional form:
   - get_fct_form_params() → create_elasticity_functions() → solve_FP()
   - Generates: Simulated_data_[spec]_alt.mat

2. GMM ESTIMATION (if run_gmm=1):
   For each specification and estimator combination:
   - Script_GMM_simulation()
   - Generates: Simulated_[spec]_Estimates_[basis]_[knots]_alt.mat

3. QQ ESTIMATION (if run_qq=1):
   For each specification:
   - Script_QQ_simulation()
   - Generates: Simulated_[spec]_Estimates_QQLogNormal_1_[force_n]_alt.mat

4. PLOTTING:
   - figures_gen(): Creates main elasticity comparison plots
   - generate_QQ_plot(): Creates QQ plot of simulated sales quantiles against the estimated ones by the QQ

All estimation and plotting functions support parallel processing and selective simulation 
running for computational efficiency.