# Data and Code for: From Heterogeneous Firms to Heterogeneous Trade Elasticities: The Aggregate Implications of Firm Export Decisions
# Replication Guide

## Adao (Chicago Booth/NBER), Arkolakis (Yale/NBER), Ganapati (Georgetown/NBER)
### Journal of Political Economy - Prepared February 2026 (updated March 2026)



## For Replicators
The code runs in two parts. One constructs code in Stata. One runs computations in Matlab. Conditional on having intermediate files, these stages can be run separately. The program MasterProgram.do is a Stata script that calls all programs in sequence. If that is not possible, then the Matlab programs (as listed below) can be run sequentially after the Stata programs are run.

## Overview
To run the Matlab programs Master_script_gravity_v3, the folder `/Data/Int/WIOD_sampleB/` must be fully populated with `.csv` for all economic fundamentals. This is done by running ALL the Stata scripts. Then to run the Matlab program `Quantiles_Colombia_Wrapper.`, additionally the files in `/Data/WorldBankData/Clean data/` must be populated by the routine `Code/p02_quantiles_colombia_firmexports.do`. See below for file path details and requirements. 
This is the second version of this guide. The first version did not anticipate case-sensitive file systems. 

## Data Availability and Provenance Statements

The author of the manuscript have legitimate access to and permission to use the data used in this manuscript. We do not have rights to redistribute all the data. Some data cannot be made publicly available. Instructions on access are detailed below.


## 1. Data Sources

The analysis relies on the following data sources. Users should obtain these files and place them in the `Data/` directory as specified.

| Data Source | Description | Location in Project | Access Information |
|-------------|-------------|---------------------|--------------------|
| **UNCTAD TRAINS** | Tariff data used for gravity estimation. | Public Use: `Data/TRAINS/` | Available via WITS (World Integrated Trade Solution) at [wits.worldbank.org](https://wits.worldbank.org/) and requires registration. We use data sourced from Boehm et al. (2023) Raw files expected in `Data/TRAINS/raw_data/` (e.g., `TRAINS1995.zip` to `TRAINS2018.zip`). |
| **World Bank EDD** | Exporter Dynamics Database. | Public Use: `Data/EDD_Data/` | Available from the World Bank at [worldbank.org/en/research/brief/exporter-dynamics-database](https://www.worldbank.org/en/research/brief/exporter-dynamics-database). Microdata `CY_manuf.dta` and `CYD_manuf.dta`. |
| **World Bank EDD - Colombian Microdata** | Colombian Official Exporter Records. Sourced Through the World Bank | Public Use: `Data/WorldBankData/` | The specific file `COL_EXP_2012.dta` is required in `Data/WorldBankData/Raw data/`. Access to this data was kindly provided by Ana Fernandes (World Bank)|
| **CEPII** | Gravity variables (distance, language, colony). | Public Use: `Data/CEPII/` | Available from CEPII (Gravity database) at [cepii.fr](http://www.cepii.fr/CEPII/en/bdd_modele/presentation.asp?id=8). Used for standard gravity controls. |
| **WIOD** | World Input-Output Database (2016 Release). | Dowload Location: `Data/WIOD/` | Available from [wiod.org](http://www.wiod.org). Used to construct the main analysis sample (`WIOD_sampleB`). |
| **Teti Tariffs** | Detailed tariff data provided by Feodora Teti. | Download Location: `Data/Teti/` | Available at [feodorateti.github.io/data.html](https://feodorateti.github.io/data.html). Files like `tariff2012_beta1-2024-12.dta` are merged in `p03_Create_Tariffs.do`. |
| **OECD SDBS/TEC** | Firm demographics and trade by enterprise characteristics. | Public Use: `Data/OECD/` | Available from OECD Stat at [stats.oecd.org](https://stats.oecd.org/). Files: `SDBS_*.csv`, `TEC3_REV4-en.csv`. Used for firm counts and exporter stats. |
| **US Census** | Number of active manufacturer enterprises in the US (2012). | Public Use: `Data/US_Census/` | Summary Series ECN_2012_US_31SG1 available at [data.census.gov](https://data.census.gov/map?q=B+B+General+Painting&tid=ECNBASIC2012.EC1231SG1&layer=VT_2012_040_00_PP_D1&loc=43.3751,-113.1138,z2.6270). |
| **Australia Data** | Characteristics of Australian Exporters. | Public Use: `Data/Australia/` | Available from the ABS at [abs.gov.au](https://www.abs.gov.au/statistics/economy/international-trade/characteristics-australian-exporters/latest-release). Includes cleaned csv (`AUS_cleaned.csv`) done by the authors. |
| **China Statistical Yearbook** | Firm statistics for China. | Public Use: `Data/China_Statistical_Yearbook/` | 2019 Yearbook available at [stats.gov.cn](https://www.stats.gov.cn/sj/ndsj/2019/indexeh.htm). Cleaned CSV made the authors.|
| **China Export Data** | Trade statistics for China. | Public Use: `Data/China_Exports/` | Export data (`Final Data.dta`) kindly provided by Peter Morrow (University of Toronto). |
| **WB Enterprise Surveys** | Firm export probabilities. | Public Use: `Data/WB_Enterprise_Data/` | Available from World Bank Enterprise Surveys at [enterprisesurveys.org](https://www.enterprisesurveys.org/). Used in `make_WB_Enterprise`. |
| **GSP Rates** | Generalized System of Preferences and MFN tariff rates (2012). | Public Use: `Data/GSP/` | Files: `mfn_gsp_2012.xlsx`. Used for counterfactual analysis. Constructed manually from the TRAINS dataset by the authors. |
| **Eora MRIO** | Multi-Region Input-Output data. | Download Location: `Data/EORA/` | Available from [worldmrio.com](https://worldmrio.com). Used to complete missing flows in WIOD. |

## 2. Programs

The replication is managed by a master Stata script (`Code/MasterProgram.do`) which calls subsequent Stata and MATLAB scripts.

### Top-Level Script
- **`Code/MasterProgram.do`**: Stata script that sets global paths and executes the cleaning, estimation, and simulation routines in order.
  - *Note:* You must set the `ROOT` global in this file to your local path before running. Additionally you must set up a path for Matlab.
  - *Note:* This assumes you have Stata 16/MP and Matlab 2023a install on *nix or OS X environment, as well as an c or z-shell terminal enivornment. We cannot provide assistance for windows users as file paths and shell scripts may not be compatible. 

### Sequence of execution (called by MasterProgram.do):

0. **Preliminary Steps**
    -   Make sure all zipped files are unzipped and located at their expected locations. In particular `/Data/TRAINS/raw_data.zip` needs to be unzipped to `/Data/TRAINS/raw_data/`. And all the files in `$ROOT/Data/Teti/tariff*_beta1-2024-12.zip`.

1.  **Data Cleaning & Construction (Stata)**
    -   `Code/p01_TRAINS_clean_data_v1_AAG.do`: Cleans raw TRAINS tariff data.
    -   `Code/p02_quantiles_colombia_firmexports.do`: Processes World Bank EDD data for Colombia firm export quantiles.
    -   `Code/p02_eora_trade_matrix_clean.do`: Cleans Eora MRIO data.
    -   `Code/p02_ChinaData_Morrow.do`: Cleans China export data.
    -   `Code/p03_Create_Tariffs.do`: Merges TRAINS data with Teti tariff data.
    -   `Code/p04_Sample_Creation.do`: Main sample construction (merges Trade, Gravity, and Production data).
    -   `Code/p05_Sample_Creation_HS2.do`: Creates HS2-level datasets for heterogeneity analysis.
    -   `Code/p06_ReducedForm.do`: Runs reduced-form gravity regressions.

2.  **Estimation & Simulation (MATLAB)**
    -   `Code/GMM_estimation_v3/Master_script_gravity_v3.m`: Main GMM estimation of gravity elasticities.
    -   `Code/GFT_v2/GFT_script_v17.m`: General Equilibrium analysis, calculating welfare gains and counterfactuals.
    -   `Code/LinearCF_v3/LinearCF_wrapper_v3.m`: Linear counterfactual simulations (Tariff Cuts, GSP removal).
    -   `Code/ColombiaTests/Quantiles_Colombia_Wrapper.m`: Monte Carlo simulations calibrated to Colombia data.
    -   `Code/MonteCarlo/run_mc_master.m`: Additional Monte Carlo baseline generation.
    -   `Code/LogCorrected/plots.m` & `Code/Appendix_LogPareto/AppendixFigure_LogPareto.m`: Robustness checks and Appendix figures for Log-Pareto specifications.

### Computational Requirements
These programs were run on Stata 16/MP and Matlab 2023a in an OS X environment (M2 processor with 16 cores and 32 GB) and a 1 TB SSD for all the programs except `run_mc_master' which was run on an external Linux cluster with 100 available cores (32 GB addressable each). Programs do require z-shell scripting, so they do need an *nix/OSX environment. All programs were tested on a case-insensitive file system. We have done our best to try to align this to run in a case-sensitive file system.

### Computational Details

Stata 16 code requires the user to install the estout, gtools, ftools, texsave, moremata, and reghdfe packages (and all their respective dependencies). Matlab requires the Statistics and Optimization toolboxes. We use the now depreciated "readcsv" and "writecsv" functions that do not work on newer Matlab versions. 

| Package | Version | Date |
|---------|---------|------|
| estout | 3.31 | 26apr2022 |
| gtools | 1.10.1 | 05Dec2022 |
| ftools | 2.49.1 | 08aug2023 |
| texsave | 1.6.1 | 22may2023 |
| moremata | (Updated with Stata 16 ) | — |
| require | (Updated with Stata 16 ) | — |
| reghdfe | 6.12.3 | 08aug2023 |



### Summary of required storage space
1TB. If space is a premium the unzipped raw tariff files (detailed below) can be deleted or zipped after the data cleaning is finished. If so, replication can be done on a tighter 250GB system. Further reduction can be done where files after p06 do not need most of the raw datasets. In this case, replication can be done on 25GB systems. However, we do not all such delete the files for completeness. 

### Summary time to reproduce
8-24 hours

### Controlled Randomness

The program `Code/MonteCarlo/run_mc_master.m` conducts Monte Carlo exercises and sets random seeds deterministically in subprograms.

## 2.5 Generated Datasets

This section documents the intermediate and final datasets generated by each script in the pipeline. The central handoff point between the Stata data-construction stage and the MATLAB estimation stage is the directory `Data/Int/WIOD_sampleB/`, which contains all CSV files used by the MATLAB programs.

### Stata: Data Cleaning and Construction

| Script | Key Datasets Generated | Used By |
|--------|----------------------|---------|
| `p01_TRAINS_clean_data_v1_AAG.do` | `Data/TRAINS/cleaned_data/T[YEAR].dta` (cleaned tariff files, 2000–2014) | `p03` |
| `p02_quantiles_colombia_firmexports.do` | `Data/WorldBankData/Clean data/colombia_exp_percentiles.csv`, `colombia_exp_names.csv` | `Quantiles_Colombia_Wrapper.m` |
| `p02_eora_trade_matrix_clean.do` | `Data/EORA/Output/trade_flows_all_sectors.csv` | `p04` |
| `p02_ChinaData_Morrow.do` | `Data/China_Exports/Final/Final Data.dta` | `p04` |
| `p03_Create_Tariffs.do` | `Data/TRAINS/collapsed/T_all_collapse.dta`, `Teti_all_collapse.dta`, `T_all_collapse_hs1.dta`, `Teti2012_all_collapse_hs1.dta` | `p04`, `p05`, `p06` |
| `p04_Sample_Creation.do` | See **WIOD_sampleB datasets** and **other intermediate datasets** tables below | `p05`, `p06`, all MATLAB scripts |
| `p05_Sample_Creation_HS2.do` | `Data/Int/WIOD_sampleB/[YEAR]_Tuw_h[H].csv`, `[YEAR]_Tetiuw_h[H].csv` (HS1 sectoral samples, H=0,3–17) | `Master_script_gravity_v3.m` (sectoral estimation) |
| `p06_ReducedForm.do` | `Output/TableOA4_ReducedFormControlsInteracted.tex` | Paper Table OA.4 |

### Datasets in `Data/Int/WIOD_sampleB/` (generated by `p04_Sample_Creation.do`)

These CSV files are the primary inputs to all MATLAB estimation scripts. Each file is indexed by year (typically `2012B`).

| File Pattern | Description | Used By |
|-------------|-------------|---------|
| `[YEAR]_Tetiuw.csv` | Baseline gravity sample (bilateral trade, tariffs from Teti, gravity controls) | `Master_script_gravity_v3.m`, `GFT_script_v17.m`, `LinearCF.m` |
| `[YEAR]_Tuw.csv` | Alternative gravity sample (tariffs from TRAINS) | `Master_script_gravity_v3.m` |
| `[YEAR]_TuwIV2.csv` | IV specification sample | `Master_script_gravity_v3.m` |
| `n_ij_[YEAR]B.csv` | Bilateral number of exporters matrix (C×C) | `GFT_script_v17.m`, `LinearCF.m` |
| `x_ij_[YEAR]B.csv` | Bilateral average export size matrix (C×C) | `GFT_script_v17.m` |
| `XX_ij_[YEAR]B.csv` | Bilateral total exports matrix (C×C) | `GFT_script_v17.m`, `LinearCF.m` |
| `balance_i_[YEAR]B.csv` | Country-level aggregates (expenditure, production, trade balance) | `GFT_script_v17.m` |
| `G_ij_[YEAR]B.csv` | Bilateral group classification (developed/developing origin) | `GFT_script_v17.m`, `LinearCF.m` |
| `G4_ij_[YEAR]B.csv` | Four-way group classification | `GFT_script_v17.m` |
| `Sample_[YEAR]B.csv` | Sample selection indicators | `GFT_script_v17.m`, `LinearCF_wrapper_v3.m` |
| `l_i_[YEAR]B.csv` | Country labels | `GFT_script_v17.m` |
| `N_ii_[YEAR]B.csv` | Domestic firm counts | `LinearCF.m` |
| `GSP_delta_[YEAR]B.csv` | GSP tariff reduction indicators | `LinearCF.m` (GSP counterfactual) |
| `GSP_[YEAR]B.csv` | GSP eligibility indicators | `LinearCF.m` |
| `mfn_rate_[YEAR]B.csv` | MFN tariff rates | `LinearCF.m` |
| `[YEAR]_survival_.csv`, `[YEAR]_survival_all_.csv`, `[YEAR]_survival_alt3_.csv` | Firm survival rates (1-year, all, 3-year) | `Master_script_gravity_v3.m` (robustness) |
| `[YEAR]_nonii_.csv` | Sample dropping imputed domestic firms | `Master_script_gravity_v3.m` (robustness) |

### Other Intermediate Datasets (generated by `p04_Sample_Creation.do`)

| File | Description | Used By |
|------|-------------|---------|
| `Data/Int/stack_data_B.dta` | Stacked panel dataset (all years, all bilateral pairs with tariffs and gravity) | `p06_ReducedForm.do` |
| `Data/Int/CEPII_Grav_temp.dta` | Cleaned CEPII gravity variables | `p05_Sample_Creation_HS2.do` |
| `Data/Int/CEPII_rich.dta` | Developed/developing country classification | `p05_Sample_Creation_HS2.do` |
| `Output/[YEAR]_Sample_summary.csv` (`.tex`) | Estimation sample summary statistics | Paper Table OA.2 |
| `Output/Hist_[YEAR]B.pdf` | Histogram of exporter firm shares | Paper Figure OA.6 |

### MATLAB: Estimation and Counterfactual Datasets

| Script | Key Datasets Generated | Used By |
|--------|----------------------|---------|
| `Master_script_gravity_v3.m` | `Code/GMM_estimation_v3/gamma_gravity.mat` (estimated gravity coefficients) | `GFT_script_v17.m`, `LinearCF.m` (both load this file) |
| `LinearCF.m` | `Code/LinearCF_v3/ESTIMATES_LINEAR.mat`, `ESTIMATES_SPLINE2.mat`, `ESTIMATES_SPLINE4.mat`, `ESTIMATES_SPLINE6.mat` (counterfactual estimation objects) | `Graphs_TC.m`, `Graphs_GSP.m` |
| `LinearCF_wrapper_v3.m` → `Graphs_TC.m` | `Output/Decomposition_Results_TC[scenario].csv` | Paper Table 1 |
| `LinearCF_wrapper_v3.m` → `Graphs_GSP.m` | `Output/Decomposition_Results_GSP[scenario].csv` | Paper Table OA.6 |
| `GFT_script_v17.m` | `Output/GFT_results_revision_2012B.xlsx` (welfare decomposition, sheets: "Generalized Pareto", "SplineOD") | Paper Figure 5 results |
| `run_mc_master.m` → `solve_FP.m` | `Code/MonteCarlo/output/data/Simulated_data_[spec]_alt.mat` (simulated economies for LogNormal, MPareto, MEstimates) | `Script_GMM_simulation.m`, `Script_QQ_simulation.m` |
| `run_mc_master.m` → `Script_GMM_simulation.m` | `Code/MonteCarlo/output/estimation/Simulated_[spec]_Estimates_[basis]_[knots]_alt.mat` | `figures_gen.m` |
| `run_mc_master.m` → `Script_QQ_simulation.m` | `Code/MonteCarlo/output/estimation/Simulated_[spec]_Estimates_QQLogNormal_1_[force_n]_alt.mat` | `figures_gen.m` |

### Dataset Dependency Chain

```
Raw Data Sources (Data/[source]/)
    │
    ├── p01 ──→ Data/TRAINS/cleaned_data/T[YEAR].dta
    │               │
    ├── p02s ─→ Data/WorldBankData/Clean data/*.csv     ──→ Quantiles_Colombia_Wrapper.m
    │         → Data/China_Exports/Final/Final Data.dta ─┐
    │         → Data/EORA/Output/*.csv                  ─┤
    │               │                                    │
    ├── p03 ──→ Data/TRAINS/collapsed/*.dta              │
    │               │                                    │
    ├── p04 ──→ Data/Int/WIOD_sampleB/*.csv  ←───────────┘
    │           Data/Int/stack_data_B.dta
    │               │
    ├── p05 ──→ Data/Int/WIOD_sampleB/*_h[H].csv (sectoral)
    │               │
    └── p06 ──→ Output/TableOA4*.tex
                    │
                    ▼
    Master_script_gravity_v3.m ──→ gamma_gravity.mat
                    │                    │
                    ├────────────────────┤
                    ▼                    ▼
    GFT_script_v17.m          LinearCF_wrapper_v3.m
         │                          │
         ▼                          ▼
    Output/GFT_results*.xlsx   Output/Decomposition_Results*.csv

    run_mc_master.m (self-contained)
         │
         ├──→ Simulated_data_[spec]_alt.mat
         ├──→ Simulated_[spec]_Estimates_*.mat
         └──→ Output figures (Elasticity_*.pdf, GFT_*.pdf)
```

## 3. Output Files

The following table maps the figures and tables in `AAG_Feb2026.pdf` to the output files and the code that generates them.

| Figure / Table | Title | Output Filename | Generating Program (Top Level -> Sub-script) |
|----------------|-------|-----------------|---------------------------------------------|
| **Figure 1** | Theoretical Trade Elasticities with Log-Normal Heterogeneity | [fig/F1_plot_theta_e_vs_log_n.eps](Draft/fig/F1_plot_theta_e_vs_log_n.eps)<br>[fig/F1_plot_theta_c_vs_log_n.eps](Draft/fig/F1_plot_theta_c_vs_log_n.eps)<br>[fig/F1_plot_theta_vs_log_n.eps](Draft/fig/F1_plot_theta_vs_log_n.eps) | `Code/MasterProgram.do` -> `Code/LogCorrected/plots.m` |
| **Figure 2** | Recovering Trade Elasticities: Monte Carlo Simulations | [fig/Elasticity_MPareto.pdf](Draft/fig/Elasticity_MPareto.pdf)<br>[fig/Elasticity_Lognormal.pdf](Draft/fig/Elasticity_Lognormal.pdf)<br>[fig/Elasticity_MEstimates.pdf](Draft/fig/Elasticity_MEstimates.pdf) | `Code/MasterProgram.do` -> `Code/MonteCarlo/run_mc_master.m` -> `Code/MonteCarlo/figures_gen.m` |
| **Figure 3** | Semiparametric Gravity of Firm Exports – Single Group | [fig/F2a_base_extensive.eps](Draft/fig/F2a_base_extensive.eps)<br>[fig/F2a_base_rho.eps](Draft/fig/F2a_base_rho.eps)<br>[fig/F2a_base_theta.eps](Draft/fig/F2a_base_theta.eps) | `Code/MasterProgram.do` -> `Code/GMM_estimation_v3/Master_script_gravity_v3.m` |
| **Figure 4** | Bilateral Trade Elasticity – Developed and Developing Origins | [fig/F3a_origin_theta_split.eps](Draft/fig/F3a_origin_theta_split.eps)<br>[fig/F3d1_Cross_theta_quadA.eps](Draft/fig/F3d1_Cross_theta_quadA.eps)<br>[fig/F3d1_Cross_theta_quadB.eps](Draft/fig/F3d1_Cross_theta_quadB.eps) | `Code/MasterProgram.do` -> `Code/GMM_estimation_v3/Master_script_gravity_v3.m` |
| **Figure 5** | The Gains From Trade: Comparison to Other Parametric Forms | [fig/relative_12345_GFT.pdf](Draft/fig/relative_12345_GFT.pdf)<br>[fig/relative_12345_GFT_LN.pdf](Draft/fig/relative_12345_GFT_LN.pdf) | `Code/MasterProgram.do` -> `Code/GFT_v2/GFT_script_v17.m` |
| **Figure 6** | Impact of a Uniform Reduction in Bilateral Trade Costs on Welfare and its Components: The Role of Parametric Assumptions | [fig/trelative_12345_1pc.pdf](Draft/fig/trelative_12345_1pc.pdf)<br>[fig/trelative_123_1pc.pdf](Draft/fig/trelative_123_1pc.pdf)<br>[fig/trelative_45_1pc.pdf](Draft/fig/trelative_45_1pc.pdf) | `Code/MasterProgram.do` -> `Code/LinearCF_v3/LinearCF_wrapper_v3.m` -> `Code/LinearCF_v3/Graphs_TC.m` |
| **Table 1** | Impact of Uniform Reductions in Trade Costs on Welfare and its Components | [fig/Decomposition_Results_TC2.csv](Draft/fig/Decomposition_Results_TC2.csv) Note: Columns 2 onward need to converted to percentages to match paper. | `Code/MasterProgram.do` -> `Code/LinearCF_v3/LinearCF_wrapper_v3.m` |
| **Figure OA.1** | Distributional Assumptions and the Firm Export Margins | [fig/AF1_pdf_e.eps](Draft/fig/AF1_pdf_e.eps)<br>[fig/AF1_fEr.eps](Draft/fig/AF1_fEr.eps) | `Code/MasterProgram.do` -> `Code/Appendix_LogPareto/AppendixFigure_LogPareto.m` |
| **Figure OA.2** | Distributional Assumptions and the Elasticity of Firm Export Margins | [fig/Extensive_v8.eps](Draft/fig/Extensive_v8.eps)<br>[Draft/fig/Thetac_v8.eps](Draft/fig/Thetac_v8.eps)<br>[fig/Theta_v8.eps](Draft/fig/Theta_v8.eps) | `Code/MasterProgram.do` -> `Code/Appendix_LogPareto/AppendixFigure_LogPareto.m` |
| **Figure OA.3** | Monte Carlo: Distributions of Gains From Trade | [fig/GFT_MPareto.pdf](Draft/fig/GFT_MPareto.pdf)<br>[fig/GFT_LogNormal.pdf](Draft/fig/GFT_LogNormal.pdf)<br>[fig/GFT_MEstimates.pdf](Draft/fig/GFT_MEstimates.pdf) | `Code/MasterProgram.do` -> `Code/GFT_v2/GFT_script_v17.m` |
| **Figure OA.4** | Monte Carlo: GMM Estimation with functional for restrictions | [fig/Elasticity_MPareto_apx.pdf](Draft/fig/Elasticity_MPareto_apx.pdf)<br>[fig/Elasticity_LogNormal_apx.pdf](Draft/fig/Elasticity_LogNormal_apx.pdf)<br>[fig/Elasticity_MEstimates_apx.pdf](Draft/fig/Elasticity_MEstimates_apx.pdf) | `Code/MasterProgram.do` -> `Code/MonteCarlo/run_mc_master.m` -> `Code/MonteCarlo/figures_gen.m` |
| **Figure OA.5** | Monte Carlo: Example QQ Estimator | [fig/QQ_scatter_MEstimates_LogNormal_k1_sim1_i1_j4.pdf](Draft/fig/QQ_scatter_MEstimates_LogNormal_k1_sim1_i1_j4.pdf) | `Code/MasterProgram.do` -> `Code/MonteCarlo/run_mc_master.m` -> `Code/MonteCarlo/figures_gen.m` |
| **Figure OA.6** | Empirical Distribution of Exporter Firm Shares, 2012 | [fig/Hist_2012B.pdf](Draft/fig/Hist_2012B.pdf) | `Code/MasterProgram.do` -> `Code/p04_Sample_Creation.do` |
| **Table OA.1** | Estimation Data Sources | No Data |
| **Table OA.2** | Estimation Data Summary | [fig/2012_Sample_summary.tex](Draft/fig/2012_Sample_summary.tex)  |  `Code/MasterProgram.do` -> `Code/p04_Sample_Creation.do`|
| **Table OA.3** | Manufacturing Sectoral Aggregates | No Data |
| **Figure OA.7** | Semiparametric Gravity of Firm Exports versus Log-Pareto Distribution of Entry Potentials | [fig/FA_logPareto_Sim.eps](Draft/fig/FA_logPareto_Sim.eps) | `Code/MasterProgram.do` -> `Code/Appendix_LogPareto/AppendixFigure_LogPareto.m` |
| **Figure OA.8** | Empirical Distribution of Bilateral Trade Elasticities in 2012 | [fig/theta_dist_all2.eps](Draft/fig/theta_dist_all2.eps) | `Code/MasterProgram.do` -> `Code/GMM_estimation_v3/Master_script_gravity_v3.m` |
| **Table OA.4** | Reduced-Form Gravity Specification |  [fig/TableOA4_ReducedFormControlsInteracted.tex](Draft/fig/TableOA4_ReducedFormControlsInteracted.tex) |  `Code/MasterProgram.do` -> `Code/p06_ReducedForm.do` |
| **Figure OA.9** | Elasticity of Firm Exports and Distributional Assumptions – Single Group | [fig/F6_lit_theta.eps](Draft/fig/F6_lit_theta.eps) | `Code/MasterProgram.do` -> `Code/GMM_estimation_v3/Master_script_gravity_v3.m` |
| **Figure OA.10** | Log-Normal QQ Estimator: Exports of Colombian Firms to the United States | [fig/QQ_Colombia_Quantiles_BaselineP.eps](Draft/fig/QQ_Colombia_Quantiles_BaselineP.eps)<br>[fig/QQ_Colombia_Quantiles_LNP.eps](Draft/fig/QQ_Colombia_Quantiles_LNP.eps) | `Code/MasterProgram.do` -> `Code/ColombiaTests/Quantiles_Colombia_Wrapper.m` |
| **Table OA.5** | Fit of QQ Estimator: Exports of Colombian Firms by Destination | [fig/QQ_Colombia_R2_Summary.csv](Draft/fig/QQ_Colombia_R2_Summary.csv) | `Code/MasterProgram.do` -> `Code/ColombiaTests/Quantiles_Colombia_Wrapper.m` |
| **Figure OA.11** | Constant-Elasticity Gravity – Developed and Developing Origins | [fig/F3a_origin_one_theta_split.eps](Draft/fig/F3a_origin_one_theta_split.eps) | `Code/MasterProgram.do` -> `Code/GMM_estimation_v3/Master_script_gravity_v3.m` |
| **Figure OA.12** | Semiparametric Gravity of Firm Exports – Origin's Income Level | [fig/F3d5_4origin_theta_quad.eps](Draft/fig/F3d5_4origin_theta_quad.eps) | `Code/MasterProgram.do` -> `Code/GMM_estimation_v3/Master_script_gravity_v3.m` |
| **Figure OA.13** | Semiparametric Gravity of Firm Exports – Within-Sector Estimation | [fig/F4_HS_all_theta.eps](Draft/fig/F4_HS_all_theta.eps)<br>[fig/F4b_HS_origin_theta_split.eps](Draft/fig/F4b_HS_origin_theta_split.eps)<br>[fig/F4b_HS_destination_theta_split.eps](Draft/fig/F4b_HS_destination_theta_split.eps) | `Code/MasterProgram.do` -> `Code/GMM_estimation_v3/Master_script_gravity_v3.m` -> `Code/GMM_estimation_v3/hs_routine_gravity.m` |
| **Figure OA.14** | Semiparametric Gravity of Firm Exports – Extensive Margin Elasticity (Sectoral) | [fig/F6_HS_3extensive.eps](Draft/fig/F6_HS_3extensive.eps) and others | `Code/MasterProgram.do` -> `Code/GMM_estimation_v3/Master_script_gravity_v3.m` -> `Code/GMM_estimation_v3/hs_routine_gravity.m` |
| **Figure OA.15** | Semiparametric Gravity of Firm Exports – Firm Composition Elasticity (Sectoral) | [fig/F6_HS_3rho.eps](Draft/fig/F6_HS_3rho.eps) and others | `Code/MasterProgram.do` -> `Code/GMM_estimation_v3/Master_script_gravity_v3.m` -> `Code/GMM_estimation_v3/hs_routine_gravity.m` |
| **Figure OA.16** | Semiparametric Gravity of Firm Exports – Bilateral Trade Elasticity (Sectoral) | [fig/F6_HS_3theta.eps](Draft/fig/F6_HS_3theta.eps) and others | `Code/MasterProgram.do` -> `Code/GMM_estimation_v3/Master_script_gravity_v3.m` -> `Code/GMM_estimation_v3/hs_routine_gravity.m` |
| **Figure OA.17** | Semiparametric Gravity of Firm Exports – Determinants of Market Integration | [fig/F3d4_Deep_theta_split.eps](Draft/fig/F3d4_Deep_theta_split.eps)<br>[fig/F3d5_LC_theta_split.eps](Draft/fig/F3d5_LC_theta_split.eps) | `Code/MasterProgram.do` -> `Code/GMM_estimation_v3/Master_script_gravity_v3.m` |
| **Figure OA.18** | Semiparametric Gravity of Firm Exports – Alternative Cost Pass-Through | [fig/F5_basekf5theta.eps](Draft/fig/F5_basekf5theta.eps)<br>[fig/F5_basekf1theta.eps](Draft/fig/F5_basekf1theta.eps)<br>[fig/F5_sigma24theta.eps](Draft/fig/F5_sigma24theta.eps)<br>[fig/F5_sigma34theta.eps](Draft/fig/F5_sigma34theta.eps) | `Code/MasterProgram.do` -> `Code/GMM_estimation_v3/Master_script_gravity_v3.m` |
| **Figure OA.19** | Semiparametric Gravity of Firm Exports, Alternative Tariff Database | [fig/F5_T2012_v_TETIextensive.eps](Draft/fig/F5_T2012_v_TETIextensive.eps)<br>[fig/F5_T2012_v_TETIrho.eps](Draft/fig/F5_T2012_v_TETIrho.eps)<br>[fig/F5_T2012_v_TETItheta.eps](Draft/fig/F5_T2012_v_TETItheta.eps) | `Code/MasterProgram.do` -> `Code/GMM_estimation_v3/Master_script_gravity_v3.m` |
| **Figure OA.20** | Semiparametric Gravity of Firm Exports – Tariff IV | [fig/F5_IV_extensive.eps](Draft/fig/F5_IV_extensive.eps)<br>[fig/F5_IV_rho.eps](Draft/fig/F5_IV_rho.eps)<br>[fig/F5_IV_theta.eps](Draft/fig/F5_IV_theta.eps) | `Code/MasterProgram.do` -> `Code/GMM_estimation_v3/Master_script_gravity_v3.m` |
| **Figure OA.21** | Semiparametric Gravity of Firm Exports – Alternative Inference | [fig/F2a_base_extensive_CC.eps](Draft/fig/F2a_base_extensive_CC.eps)<br>[fig/F2a_base_rho_CC.eps](Draft/fig/F2a_base_rho_CC.eps) | `Code/MasterProgram.do` -> `Code/GMM_estimation_v3/Master_script_gravity_v3.m` -> `Code/GMM_estimation_v3/CC_se.m` |
| **Figure OA.22** | Semiparametric Gravity of Firm Exports – Alternative Functional Form | [fig/F4_k4_extensive.eps](Draft/fig/F4_k4_extensive.eps)<br>[fig/F4_k4_rho.eps](Draft/fig/F4_k4_rho.eps)<br>[fig/F4_k4_theta.eps](Draft/fig/F4_k4_theta.eps) | `Code/MasterProgram.do` -> `Code/GMM_estimation_v3/Master_script_gravity_v3.m` |
| **Figure OA.23** | Semiparametric Gravity of Firm Exports – Single Group, 2010 | [fig/F5_y2010extensive.eps](Draft/fig/F5_y2010extensive.eps)<br>[fig/F5_y2010rho.eps](Draft/fig/F5_y2010rho.eps)<br>[fig/F5_y2010theta.eps](Draft/fig/F5_y2010theta.eps) | `Code/MasterProgram.do` -> `Code/GMM_estimation_v3/Master_script_gravity_v3.m` |
| **Figure OA.24** | Semiparametric Gravity of Firm Exports – Single Group, 2014 | [fig/F5_y2014extensive.eps](Draft/fig/F5_y2014extensive.eps)<br>[fig/F5_y2014rho.eps](Draft/fig/F5_y2014rho.eps)<br>[fig/F5_y2014theta.eps](Draft/fig/F5_y2014theta.eps) | `Code/MasterProgram.do` -> `Code/GMM_estimation_v3/Master_script_gravity_v3.m` |
| **Figure OA.25** | Semiparametric Gravity of Firm Exports – 3-year Survival Rate | [fig/F5_nii3extensive.eps](Draft/fig/F5_nii3extensive.eps)<br>[fig/F5_nii3rho.eps](Draft/fig/F5_nii3rho.eps)<br>[fig/F5_nii3theta.eps](Draft/fig/F5_nii3theta.eps) | `Code/MasterProgram.do` -> `Code/GMM_estimation_v3/Master_script_gravity_v3.m` |
| **Figure OA.26** | Semiparametric Gravity of Firm Exports – Dropping Observations with Imputed $n_{ii}$ | [fig/F5_noimpniiextensive.eps](Draft/fig/F5_noimpniiextensive.eps)<br>[fig/F5_noimpniirho.eps](Draft/fig/F5_noimpniirho.eps)<br>[fig/F5_noimpniitheta.eps](Draft/fig/F5_noimpniitheta.eps) | `Code/MasterProgram.do` -> `Code/GMM_estimation_v3/Master_script_gravity_v3.m` |
| **Figure OA.27** | Gains from Trade: Entry and Selection of Domestic Firms | [fig/2_groups_n_Ni.pdf](Draft/fig/2_groups_n_Ni.pdf)<br>[fig/2_groups_n_nii.pdf](Draft/fig/2_groups_n_nii.pdf) | `Code/MasterProgram.do` -> `Code/GFT_v2/GFT_script_v17.m` |
| **Figure OA.28** | Gains from Trade: Firm Profit Margins | [fig/2_groups_p2_nii.pdf](Draft/fig/2_groups_p2_nii.pdf)<br>[fig/2_groups_p_nii.pdf](Draft/fig/2_groups_p_nii.pdf) | `Code/MasterProgram.do` -> `Code/GFT_v2/GFT_script_v17.m` |
| **Figure OA.29** | Impact of a Uniform Reduction in Trade Costs on Welfare and its Components | [fig/Decomposition1pc.pdf](Draft/fig/Decomposition1pc.pdf) | `Code/MasterProgram.do` -> `Code/LinearCF_v3/LinearCF_wrapper_v3.m` |
| **Figure OA.30** | Impact of a Uniform Reduction in Trade Costs on Firm Entry and Selection | [fig/Ni_Hat1pc.pdf](Draft/fig/Ni_Hat1pc.pdf)<br>[fig/nji2_Hat1pc.pdf](Draft/fig/nji2_Hat1pc.pdf)<br>[fig/nji_Hat1pc.pdf](Draft/fig/nji_Hat1pc.pdf)<br>[fig/nii_Hat1pc.pdf](Draft/fig/nii_Hat1pc.pdf) | `Code/MasterProgram.do` -> `Code/LinearCF_v3/LinearCF_wrapper_v3.m` |
| **Figure OA.31** | Impact of a Uniform Reduction in Trade Costs on Firm Export Margins: The Role of Parametric Assumptions | [fig/n1pc.pdf](Draft/fig/n1pc.pdf)<br>[fig/xbar_1pc.pdf](Draft/fig/xbar_1pc.pdf)<br>[fig/X_1pc.pdf](Draft/fig/X_1pc.pdf) | `Code/MasterProgram.do` -> `Code/LinearCF_v3/LinearCF_wrapper_v3.m` |
| **Table OA.6** | Impact of Heterogeneous Reductions in Trade Costs on Welfare and its Components (GSP) | [fig/Decomposition_Results_GSP*.csv](fig/) | `Code/MasterProgram.do` -> `Code/LinearCF_v3/LinearCF_wrapper_v3.m` |
| **Figure OA.32** | Impact of Reducing the Cost of Exporting from Developing to Developed Countries on Welfare and its Components | [fig/relative_12345_gsp_like.pdf](Draft/fig/relative_12345_gsp_like.pdf)<br>[fig/relative_123_gsp_like.pdf](Draft/fig/relative_123_gsp_like.pdf)<br>[fig/relative_45_gsp_like.pdf](Draft/fig/relative_45_gsp_like.pdf) | `Code/MasterProgram.do` -> `Code/LinearCF_v3/LinearCF_wrapper_v3.m` |
| **Figure OA.33** | Impact of Reducing the Cost of Exporting from Developing to Developed Countries on Firm Export Margins | [fig/nGSP_like.pdf](Draft/fig/nGSP_like.pdf)<br>[fig/xbar_GSP_like.pdf](Draft/fig/xbar_GSP_like.pdf)<br>[fig/X_GSP_like.pdf](Draft/fig/X_GSP_like.pdf) | `Code/MasterProgram.do` -> `Code/LinearCF_v3/LinearCF_wrapper_v3.m` |

## 4. References

Australian Bureau of Statistics. 2019. "Characteristics of Australian Exporters." https://www.abs.gov.au/statistics/economy/international-trade/characteristics-australian-exporters/latest-release (accessed April 12, 2019).

CEPII. 2026. "Gravity Database." http://www.cepii.fr/CEPII/en/bdd_modele/presentation.asp?id=8 (accessed February 1, 2017).

Eora. 2024. "Multi-Region Input-Output Database." https://worldmrio.com (accessed January 23, 2024).

National Bureau of Statistics of China. 2019. "China Export Data." Export data (`Final Data.dta`) kindly provided by Peter Morrow (University of Toronto). (accessed August 27, 2020).

National Bureau of Statistics of China. 2019. "China Statistical Yearbook (2019)." https://www.stats.gov.cn/sj/ndsj/2019/indexeh.htm (accessed September 2, 2020).

OECD. 2026. "Structural Business Statistics and Trade by Enterprise Characteristics." https://stats.oecd.org/ (accessed November 25, 2017).

Teti, Feodora. 2026. "Detailed Tariff Data." https://feodorateti.github.io/data.html (accessed July 18, 2025).

U.S. Census Bureau. 2026. "Summary Series ECN_2012_US_31SG1." https://data.census.gov/map?q=B+B+General+Painting&tid=ECNBASIC2012.EC1231SG1&layer=VT_2012_040_00_PP_D1&loc=43.3751,-113.1138,z2.6270 (accessed April 11, 2019).

UNCTAD. 2024. "Generalized System of Preferences (GSP) and MFN Rates." https://unctad.org/topic/trade-agreements/generalized-system-of-preferences (accessed September 23, 2024).

UNCTAD. 2026. "UNCTAD TRAINS Tariff Data." https://wits.worldbank.org/ (accessed October 13, 2022).

WIOD. 2026. "World Input-Output Database (2016 Release)." http://www.wiod.org (accessed September 2, 2020).

World Bank. 2019. "Enterprise Surveys." https://www.enterprisesurveys.org/ (accessed April 10, 2019).

World Bank. 2026. "Exporter Dynamics Database." https://www.worldbank.org/en/research/brief/exporter-dynamics-database (accessed December 9, 2015).

World Bank. 2026. "Exporter Dynamics Database - Colombian Microdata." https://www.worldbank.org/en/research/brief/exporter-dynamics-database (accessed December 9, 2015).
