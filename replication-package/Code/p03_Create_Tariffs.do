*********************************************************************************
*Code to create Data for estimation using WIOD Sample
*Uses Data generously from Boehm et al
*********************************************************************************

clear all
set type double, perm


/////////////////////////////
///Final dataset Collapse ///
/////////////////////////////
forvalues year=2002(2)2014 {
	use "$ROOT/Data/TRAINS/cleaned_data/T`year'.dta", clear
	cap mkdir "$ROOT/Data/TRAINS/collapsed"
	preserve
	gcollapse (mean) simpleAHS_w = simpleAHS simpleMFN_w = simpleMFN [aw=imports_trains], by(exporter importer year)
	tempfile t1
	saveold `t1', replace
	restore
	gcollapse (mean) simpleAHS_uw = simpleAHS simpleMFN_uw = simpleMFN (count) simpleAHS_obs = simpleAHS (sum) imports_trains, by(exporter importer year)
	merge 1:1 exporter importer year using `t1', nogen
	replace simpleAHS_w = simpleAHS_uw if simpleAHS_w == .
	saveold "$ROOT/Data/TRAINS/collapsed/T`year'.dta", replace
}

// Stack Tariff Data
clear
forvalues year=2000(2)2014 {
	cap append using "$ROOT/Data/TRAINS/collapsed/T`year'.dta"
}
tab year
rename (exporter importer) (iso3_o iso3_d)
label var simpleAHS_w "Tariff - Weighted"
label var simpleAHS_uw "Tariff - Unweighted"
label var simpleAHS_obs "Tariff - Observations"
saveold "$ROOT/Data/TRAINS/collapsed/T_all_collapse.dta", replace

// Merge in Teti Data
forvalues year=2010(2)2014 {
	use "$ROOT/Data/TRAINS/cleaned_data/T`year'.dta", clear
	rename exporter iso3_o
	rename importer is3_d
	rename hs6 hs92
	destring hs92, force replace
	merge 1:1 year hs92 iso3_o is3_d using "$ROOT/Data/Teti/tariff`year'_beta1-2024-12.dta", keep(master match)
	keep if _merge == 3
	rename is3_d iso3_d
	saveold "$ROOT/Data/TRAINS/collapsed/Teti_`year'_merged.dta", replace
  cap erase  "$ROOT/Data/Teti/tariff`year'_beta1-2024-12.dta"
}

// Stack Teti Data
clear
forvalues year=2010(2)2014 {
	append using "$ROOT/Data/TRAINS/collapsed/Teti_`year'_merged.dta"
}
tab year
gcollapse (mean) TsimpleAHS_uw = tariff TsimpleMFN_uw = mfn (count) TsimpleAHS_obs = tariff, by(iso3_o iso3_d year)
label var TsimpleAHS_uw "Teti Tariff - Unweighted"
label var TsimpleMFN_uw "Teti MFN Tariff - Unweighted"
label var TsimpleAHS_obs "Teti Tariff - Observations"
saveold "$ROOT/Data/TRAINS/collapsed/Teti_all_collapse.dta", replace

/* 
    01-05  Animal & Animal Products
    06-15  Vegetable Products
    16-24  Foodstuffs
    25-27  Mineral Products
    28-38  Chemicals & Allied Industries
    39-40  Plastics / Rubbers
    41-43  Raw Hides, Skins, Leather, & Furs
    44-49  Wood & Wood Products
    50-63  Textiles
    64-67  Footwear / Headgear
    68-71  Stone / Glass
    72-83  Metals
    84-85  Machinery / Electrical
    86-89  Transportation
    90-97  Miscellaneous */

// Baseline HS1
	forvalues year=2012/2012 {
		use "$ROOT/Data/TRAINS/cleaned_data/T`year'.dta", clear
		gen hs2 = substr(hs6,1,2)
		destring hs2, replace
        gen hs1 = .
        replace hs1 = 1 if inrange(hs2,01,05) //  Animal & Animal Products
        replace hs1 = 2 if inrange(hs2,06,15) //  Vegetable Products
        replace hs1 = 3 if inrange(hs2,16,24) //  Foodstuffs
        replace hs1 = 4 if inrange(hs2,25,27) //  Mineral Products
        replace hs1 = 5 if inrange(hs2,28,38) //  Chemicals & Allied Industries
        replace hs1 = 6 if inrange(hs2,39,40) //  Plastics / Rubbers
        replace hs1 = 9 if inrange(hs2,41,43) //  Raw Hides, Skins, Leather, & Furs
        replace hs1 = 8 if inrange(hs2,44,49) //  Wood & Wood Products
        replace hs1 = 9 if inrange(hs2,50,63) //  Textiles
        replace hs1 = 9 if inrange(hs2,64,67) //  Footwear / Headgear
        replace hs1 = 4 if inrange(hs2,68,71) //  Stone / Glass
        replace hs1 = 12 if inrange(hs2,72,83) //  Metals
        replace hs1 = 13 if inrange(hs2,84,84) //  Machinery / Electrical
        replace hs1 = 17 if inrange(hs2,85,85) //  Machinery / Electrical
        replace hs1 = 14 if inrange(hs2,86,89) //  Transportation
        replace hs1 = 15 if inrange(hs2,90,92) //  Scientific Equipment (This is a grab bag, hard to do right - we should drop.)
        replace hs1 = 16 if inrange(hs2,93,97) //  Miscellaneous 
        preserve
          gcollapse (mean) simpleAHS_w = simpleAHS simpleMFN_w = simpleMFN [aw=imports_trains], by(exporter importer year hs1)
          tempfile t1
          saveold `t1', replace
        restore
        gcollapse (mean) simpleAHS_uw = simpleAHS simpleMFN_uw = simpleMFN (count) simpleAHS_obs = simpleAHS (sum) imports_trains, by(exporter importer year hs1)
        merge 1:1 exporter importer year hs1 using `t1', nogen
        replace simpleAHS_w = simpleAHS_uw if simpleAHS_w == .
        saveold "$ROOT/Data/TRAINS/collapsed/T`year'_hs1.dta", replace
	}

  // Stack Teti Tariff Data at HS1
	clear
	forvalues year=2012/2012 {
		append using "$ROOT/Data/TRAINS/collapsed/T`year'_hs1.dta"
	}
	tab year
	rename (exporter importer) (iso3_o iso3_d)
	label var simpleAHS_w "Tariff - Weighted"
	label var simpleAHS_uw "Tariff - Unweighted"
	label var simpleAHS_obs "Tariff - Observations"
	saveold "$ROOT/Data/TRAINS/collapsed/T_all_collapse_hs1.dta", replace


// Teti HS1
    use "$ROOT/Data/TRAINS/collapsed/Teti_2012_merged.dta", clear
    gen hs2 = floor(hs92/1000)
    gen hs1 = .
    replace hs1 = 1 if inrange(hs2,01,05) //  Animal & Animal Products
    replace hs1 = 2 if inrange(hs2,06,15) //  Vegetable Products
    replace hs1 = 3 if inrange(hs2,16,24) //  Foodstuffs
    replace hs1 = 4 if inrange(hs2,25,27) //  Mineral Products
    replace hs1 = 5 if inrange(hs2,28,38) //  Chemicals & Allied Industries
    replace hs1 = 6 if inrange(hs2,39,40) //  Plastics / Rubbers
    replace hs1 = 9 if inrange(hs2,41,43) //  Raw Hides, Skins, Leather, & Furs
    replace hs1 = 8 if inrange(hs2,44,49) //  Wood & Wood Products
    replace hs1 = 9 if inrange(hs2,50,63) //  Textiles
    replace hs1 = 9 if inrange(hs2,64,67) //  Footwear / Headgear
    replace hs1 = 4 if inrange(hs2,68,71) //  Stone / Glass
    replace hs1 = 12 if inrange(hs2,72,83) //  Metals
    replace hs1 = 13 if inrange(hs2,84,84) //  Machinery / Electrical
    replace hs1 = 17 if inrange(hs2,85,85) //  Machinery / Electrical
    replace hs1 = 14 if inrange(hs2,86,89) //  Transportation
    replace hs1 = 15 if inrange(hs2,90,92) //  Scientific Equipment (This is a grab bag, hard to do right - we should drop.)
    replace hs1 = 16 if inrange(hs2,93,97) //  Miscellaneous 
    gcollapse (mean) simpleAHS_uw = simpleAHS simpleMFN_uw = simpleMFN TsimpleAHS_uw = tariff TsimpleMFN_uw = mfn (count) simpleAHS_obs = simpleAHS TsimpleAHS_obs = tariff (sum) imports_trains, by(iso3_o iso3_d year hs1)

    regress TsimpleAHS_uw simpleAHS_uw
    drop simpleAHS_uw simpleMFN_uw simpleAHS_obs imports_trains
    saveold "$ROOT/Data/TRAINS/collapsed/Teti2012_all_collapse_hs1.dta", replace
