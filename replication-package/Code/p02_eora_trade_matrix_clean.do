global eora "$ROOT/Data/EORA/"
global raw_data "$eora/Raw_data"
global clean_data "$eora/Clean_data"
global output "$eora/Output"

cap mkdir "$eora/Clean_data"
cap mkdir "$eora/Output"



*** All sectors

* Prepare row and column indexes (they are the same in the trade matrix)
import delimited "$raw_data/Eora26_2012_bp/labels_T.txt", clear
keep v2 v4
rename (v2 v4) (iso3_o sector)
duplicates drop
gen index = _n
tempfile labels
save `labels'

* Prepare the trade matrix in long format
import delimited "$raw_data/Eora26_2012_bp/Eora26_2012_bp_T.txt", clear
gen index = _n
merge 1:1 index using `labels', assert(match) nogen
drop index
collapse (sum) v*, by(iso3_o)
reshape long v, i(iso3_o) j(index)
rename (iso3_o v) (origin X)
merge m:1 index using `labels', assert(match) nogen
rename iso3_o destination
save "$clean_data/eora_long_all_sectors_all_cty.dta", replace

* Distinguish between RoW and countries to be included
use "$clean_data/eora_long_all_sectors_all_cty.dta", clear
drop if origin == "ROW" | destination == "ROW"
collapse (sum) X, by(origin destination)
preserve
	import delimited "$raw_data/2012_CF_sample.csv", clear
	keep if include == "to include"
	rename iso3_o origin
	tempfile origin_included
	save `origin_included'
	rename origin destination
	tempfile destination_included
	save `destination_included'
restore
merge m:1 origin using `origin_included', assert(master match) nogen
replace origin = "ROW" if include == ""
drop include
merge m:1 destination using `destination_included', assert(master match) nogen
replace destination = "ROW" if include == ""
collapse (sum) X, by(origin destination)
gen year = 2012
export delimited "$output/trade_flows_all_sectors.csv", replace


