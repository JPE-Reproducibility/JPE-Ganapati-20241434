// Export Colombia firm export percentiles
clear all

global raw_data "$ROOT/data//WorldBankData/Raw data"
global clean_data "$ROOT/data//WorldBankData/Clean data"

// Read in data from WOrld Bank
use "$raw_data/COL_EXP_2012.dta", clear
collapse (sum) v, by(d f)

gen lv = ln(v)

foreach p of numlist 1/99 {
	egen exp_p`p' = pctile (lv), by(d) p(`p')
}

// Collapse to get means, sums and counts
collapse (mean) exp_p* exp_avg=v (sum) exp_agg=v (count) exp_n=v, by(d) fast
egen tot = sum(exp_agg)
gen exp_sh = exp_agg/tot
egen exp_rank = rank(exp_agg), field 
drop tot exp_agg
sort exp_rank 
order d exp_rank exp_sh exp_n exp_avg exp_p*
rename exp_* * 


// Export data
save "$clean_data/colombia_exp_percentiles.dta", replace

preserve
	keep d
	outsheet using "$clean_data/colombia_exp_names.csv", replace
restore

preserve
	keep rank sh n avg p*
	outsheet using "$clean_data/colombia_exp_percentiles.csv", replace comma
restore