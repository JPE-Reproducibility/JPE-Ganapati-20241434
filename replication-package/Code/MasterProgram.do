// Master Program for Replication
// From Heterogeneous Firms to Heterogeneous Trade Elasticities: The Aggregate Implications of Firm Export Decisions
// See replication.md for full details on how to run the code and replicate the results in the paper

if  "`c(username)'" == "sharat" {
  global ROOT "/Users/sharat/Dropbox/Economics/2018_NonParametric_MC/Replication/"
}
if  "`c(username)'" == "user" {
  global ROOT "XXXXXX"
}

clear all
set type double, perm
global DATA   $ROOT/Data/
global INT    $ROOT/Data/Int/
global OUT    $ROOT/Output/
cap mkdir $DATA
cap mkdir $INT
cap mkdir $OUT
cap mkdir "$INT/WIOD_sampleB/"
cap mkdir "$INT/DIST/"
cap mkdir "$OUT/data/"
cap mkdir "$OUT/estimation/"
cap mkdir "$OUT/plots/"
cap mkdir "$OUT/plots/diagnostic/"
cap mkdir "$OUT/plots/additional specs/"
cap mkdir "$DATA/TRAINS/cleaned_data/"
cap mkdir "$DATA/TRAINS/collapsed/"
cap mkdir "$DATA/WorldBankData/Clean data/"


// Install packages
// We use the following packages in our code: estout, gtools, ftools, texsave, moremata, and reghdfe
// If you do not have these packages installed, please run the following code to install them before running the rest of the code
// If you have these packages installed, you can skip this step
// Note that some of these packages have dependencies on other packages, so you may need to install those as well
foreach pkg in estout require gtools ftools texsave moremata reghdfe {
    cap which `pkg'
    if _rc {
        ssc install `pkg', replace
    }
    else {
        which `pkg'
    }
}

// Check Matlab Pathing
tempfile mpath
! which matlab > `mpath' 2>&1
tempname fh
file open `fh' using `mpath', read text
file read `fh' line
file close `fh'
if strpos("`line'", "matlab") == 0 {
    di as error "MATLAB not found on PATH. Please add it to your PATH or create a symlink, e.g.:"
    di as error "  sudo ln -s /Applications/MATLAB_R2023a.app/bin/matlab /usr/local/bin/matlab"
    exit 198
}
else {
    di as result "MATLAB found at: `line'"
}

// Run Cleaning Code
// This preliminary code cleans all raw data files and creates intermediate data files for analysis
// It only needs to be run once unless raw data files are changed
// We provide all the raw files for initial replication, but do not own the rights to distribute some of them
// Therefore, some of the raw data files need to be obtained separately by the user

// Files with rights to distribute
// - TRAINS (originally from comtrade.un.org, use replication package of Boehm et al. XXX)
// - GSP data (available from XXXX)
// - Australia Firm Export Data 
// - US Firm Data 
// - China Firm Export Data (sourced via Morrow et al. XXX)
// - China Statistical Yearbook Data 
// - OECD TEC Data (sourced via OECD.stat)
// - BACI/CEPII Trade Data (sourced via www.cepii.fr)
// - World Bank EDD
// - World Bank Enterprise Surveys
// - Colombia Firm Export Data (sourced via the World Bank)

// Files without rights to distribute
// - WIOD Input-Output Tables (sourced via www.wiod.org)
// - EORA Input-Output Tables (sourced via EORA database, www.eora.global)
// - Teti Tariff Data (sourced via Teti et al. XXX)

    // TRAINS CODE 
    do $ROOT/Code/p01_TRAINS_clean_data_v1_AAG.do

    // Columbia Data
    do $ROOT/Code/p02_quantiles_colombia_firmexports.do

    // EORA Data
    do $ROOT/Code/p02_eora_trade_matrix_clean.do

    // China Data
    do $ROOT/Code/p02_ChinaData_Morrow.do

    // Run Teti Code
    do $ROOT/Code/p03_Create_Tariffs.do

// Aggregate Data Aggregation and Sample Creation
    do $ROOT/Code/p04_Sample_Creation.do
    do $ROOT/Code/p05_Sample_Creation_HS2.do

// Run Stata Analysis
    do $ROOT/Code/p06_ReducedForm.do

// Run Matlab Code

// Run GMM Wrapper
    cd $ROOT/Code/GMM_estimation_v3/
    ! matlab -nodisplay -nosplash -nodesktop -r "run('Master_script_gravity_v3.m'); exit" && echo "completed"

// Run GFT Wrapper
    cd $ROOT/Code/GFT_v2/
    ! matlab -nodisplay -nosplash -nodesktop -r "run('GFT_script_v17.m'); exit" && echo "completed"

// Run Counterfactual Wrapper
    cd $ROOT/Code/LinearCF_v3/
    ! matlab -nodisplay -nosplash -nodesktop -r "run('LinearCF_wrapper_v3.m'); exit" && echo "completed"

// Count Monte Carlo Simulations
    // Colombia Tests
    cd $ROOT/Code/ColombiaTests/
    ! matlab -nodisplay -nosplash -nodesktop -r "run('Quantiles_Colombia_Wrapper.m'); exit" && echo "completed"

// Run More Monte Carlo Simulation - Generate Baseline Economies
    cd $ROOT/Code/MonteCarlo/
    ! matlab -nodisplay -nosplash -nodesktop -r "run('run_mc_master.m'); exit" && echo "completed"

// Run Log-Corrected Simulation
    cd $ROOT/Code/LogCorrected
    ! matlab -nodisplay -nosplash -nodesktop -r "run('plots.m'); exit" && echo "completed"

// Run Log-Corrected Spline
    cd $ROOT/Code/Appendix_LogPareto
    ! matlab -nodisplay -nosplash -nodesktop -r "run('AppendixFigure_LogPareto.m'); exit" && echo "completed"



