*********************************************************************************
* Uses CEPII Data
*********************************************************************************

use  "$ROOT/Data/Int/stack_data_B", clear
egen group = group(o d)

// Just for a year
keep if year ==  2012

keep if n_ij != .
egen rank_n_ij = rank(n_ij)
replace rank_n_ij = rank_n_ij/_N
cap drop cut_n_ij
egen cut_n_ij = cut(rank_n_ij), at(0, .50, 1.01)   icodes

gen tariff_use = TsimpleAHS_uw/100
gen ltariff_use = log(1+tariff_use)

gen neglX_ij = -lX_ij
gen dist_med = cond(ldistw> 8.562,1,0)

// What We Use In Paper
global FE1 o cut_n_ij d
global FE2 o cut_n_ij#rich_orig d

global controls2  i.comcur i.contig i.colony // i.comlang_off
eststo A1: reghdfe neglX_ij c.(ldistw)##i.cut_n_ij  , absorb($FE1) cluster(group) 
estadd local FE "i,j I(n_ij)"

eststo A1ad: reghdfe neglX_ij c.(ldistw)##i.cut_n_ij $controls2, absorb($FE1) cluster(group) 
estadd local controls "X"
estadd local FE "i,j I(n_ij)"

// Double Interacted
preserve 
  cap drop rank_n_ij
  keep if n_ij != .
  bys rich_orig: egen rank_n_ij = rank(n_ij)
  bys rich_orig: egen NN = count(n_ij)

  replace rank_n_ij = rank_n_ij/NN
  cap drop cut_n_ij
  egen cut_n_ij = cut(rank_n_ij), at(0, .50, 1.01)   icodes
  eststo A4: reghdfe neglX_ij c.(ldistw )##i.cut_n_ij#i.rich_orig  , absorb($FE2) cluster(group) 
  estadd local FE "i,j I(n_ij)xDeveloped"

  eststo A4a: reghdfe neglX_ij c.(ldistw )##i.cut_n_ij#i.rich_orig  $controls2, absorb($FE2) cluster(group) 
  estadd local controls "X"
  estadd local FE "i,j I(n_ij)xDeveloped"
restore

// Output Table
global G3 se label stats( N r2 FE controls) keep(ldistw 1.cut_n_ij#c.ldistw 0.rich_orig#c.ldistw 1.cut_n_ij#0.rich_orig#c.ldistw 1.rich_orig#c.ldistw 1.cut_n_ij#1.rich_orig#c.ldistw   ) order(ldistw 1.cut_n_ij#c.ldistw 0.rich_orig#c.ldistw 1.cut_n_ij#0.rich_orig#c.ldistw 1.rich_orig#c.ldistw 1.cut_n_ij#1.rich_orig#c.ldistw  ) mtitles( "Baseline" "+ Controls" "Split" "+ Controls") replace nostar
esttab A1 A1ad A4 A4a , $G3
esttab A1 A1ad A4 A4a  using "$OUT/TableOA4_ReducedFormControlsInteracted.tex", $G3


