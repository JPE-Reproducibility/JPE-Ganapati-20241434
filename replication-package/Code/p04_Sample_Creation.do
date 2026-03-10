*********************************************************************************
*Code to create Data for estimation using Full Sample
*********************************************************************************


capture program drop MAIN
program define MAIN
  
  di "MAIN"

  FX
  make_CEPII

  setup_WIOD  
  forvalues y = 2000(1)2014 {
    di "`y'"
    global year = `y'
    make_WIOD_small, year(${year})
    make_BACI_small, year(${year})
    make_OECD_SDBS, year(${year})
    make_WB_Enterprise, year(${year})
    make_WB_EDD, year(${year})
    make_OECD_TEC, year(${year})
    make_OECD_TEC_ALL_FIRM, year(${year})
    make_AUS_exporters, year(${year})
    make_CHN_exporters, year(${year})
    make_OECD_entry, year(${year})
    stack_data, year(${year})
  } 
  

  STACK_ALL_YEARS

  global year = 2012

  GSP_2012_data
  DATA_EXPORT_AGG, year(2012)
  DATA_EXPORT_AGG, year(2014)
  DATA_EXPORT_AGG, year(2010)
  DATA_EXPORT_CF_BIG, year(2012)

end

capture program drop FX
program define FX
  import excel "$ROOT/Data/WIOD/EXR_WIOD_Nov16.xlsx", sheet("EXR") cellrange(A4:Q47) firstrow clear
  rename Acronym iso3_o
  drop Country
  rename _*2* US_FX*2*
  reshape long US_FX, i(iso) j(year)
  saveold "$ROOT/Data/Int/WIOD_FX", replace

end


// Helper Programs
capture program drop base_year_normalize
program define base_year_normalize
syntax varlist, year(integer) fe(string) prefix(string)
  di "base year: `year'"
  di "fe       : `fe'"
  di "prefix   : `prefix'"
  foreach v of varlist `varlist' {
    di "`v'"
    cap drop `prefix'`v'
    cap drop _temp*
    gen _temp = `v' if year == `year'
    bys `fe': egen _temp2 = max(_temp)
    gen `prefix'`v' = `v'/_temp2
    local label : variable label `v'
    label variable `prefix'`v' `"Normalized: `label'"'
    drop _temp _temp2
  }
  sort `fe' year
end

// Step 1
// Bilateral Trade Flows
// WIOD
capture program drop make_WIOD_small
program define make_WIOD_small
syntax, year(integer)

  use "$ROOT/Data/WIOD/wiot_stata_Nov16/WIOT`year'_October16_ROW.dta", clear
  // gen Icode1 = substr(IndustryCode,1,1)
  preserve
    keep IndustryCode IndustryDescription
    duplicates drop
    outsheet using "$ROOT/Data/Int/WIOD_sectors.csv", replace
  restore

  fcollapse (sum) v* ,by(Country) fast


  levelsof Country, local(c) clean
  foreach i in `c' "ROW" {
    di "`i'"
    * drop sum_`i'
    cap egen sum_`i' = rowtotal(v`i'*)
    cap drop v`i'*
  }
  cap drop v*
  reshape long sum_, i(Country) j(Destination) string

  rename Country iso3_o
  rename Destination iso3_d
  format  iso* %9s
  rename sum_ WIOD_X_ij
  saveold "$ROOT/Data/Int/WIOT`year'_small", replace
  keep iso3_o
  duplicates drop
  saveold "$ROOT/Data/Int/country_list_small", replace

end


capture program drop setup_WIOD 
program define setup_WIOD   



  insheet using "$ROOT/Data/CEPII/BACI_HS92_V202201/country_codes_V202201.csv",clear
  keep country_code iso_3digit_alpha
  replace iso_3digit_alpha = "TWN" if country_code == 490
  rename (country_code iso_3digit_alpha) (iso3num_o iso3_o)
  gduplicates drop
  duplicates list iso3num_o
  saveold "$ROOT/Data/Int/iso3num_o", replace

  insheet using "$ROOT/Data/CEPII/BACI_HS92_V202201/country_codes_V202201.csv",clear
  keep country_code iso_3digit_alpha
  replace iso_3digit_alpha = "TWN" if country_code == 490
  rename (country_code iso_3digit_alpha) (iso3num_d iso3_d)
  gduplicates drop
  duplicates list iso3num_d
  saveold "$ROOT/Data/Int/iso3num_d", replace



end


// STEP 1B
capture program drop make_BACI_small
program define make_BACI_small
syntax, year(integer)
  insheet using  "$ROOT/Data/CEPII/BACI_HS92_V202201/BACI_HS92_Y`year'_V202201.csv", clear
  rename (i j t) (iso3num_o iso3num_d year)
  cap destring q, force replace
  gcollapse (sum) BACI_v = v BACI_q = q, by(iso3num_o iso3num_d year)
  fmerge m:1 iso3num_o   using "$ROOT/Data/Int/iso3num_o", keep(master matched)
  assert _merge == 3
  drop _merge
  fmerge m:1 iso3num_d   using "$ROOT/Data/Int/iso3num_d", keep(master matched)
  assert _merge == 3
  drop _merge
  replace BACI_v = BACI_v/1000
  replace BACI_q = BACI_q/1000
  saveold "$ROOT/Data/Int/BACI`year'_small", replace
end

// Step 2
// Number of Active Firms
capture program drop make_OECD_SDBS
program define make_OECD_SDBS
syntax, year(integer)

    insheet using "$ROOT/Data/OECD/SDBS Business Demography Indicators/SSIS_BSC_ISIC4_11042019175351230.csv", clear
    keep if size == "Total"
    keep if var == "ENTR"
    keep if isic4 == "10_33"
    compress
    tab time
    gen year_diff = abs(time - `year')
    gsort country year_diff -time
    by country: keep if _n == 1
    // DIE
    // keep if time == `year'
    keep location   v12    value
    rename v12 SSIS_year
    rename location iso3_o
    rename value SSIS_N_ii
    cap duplicates drop
    cap saveold "$ROOT/Data/Int/SSIS`year'_small", replace

    insheet using "$ROOT/Data/OECD/SDBS Business Demography Indicators/SDBS_BDI_ISIC4_11042019204932050.csv", clear
    keep if scl == "TOTAL"
    keep if unitcode == "NBR"
    keep if ind == "ENTR_BD_EMPL"
    keep if sec == "10_33"
    gen year_diff = abs(time - `year')
    gsort country year_diff -time
    by country: keep if _n == 1
    // keep if time == `year'
    compress
    keep location country variable time value flags
    rename time SDBS_year
    rename location iso3_o
    drop flags
    duplicates drop
    gen source = "OECD SDBS"
    saveold "$ROOT/Data/Int/SDBS`year'", replace

    insheet using "$ROOT/Data/US_Census/2012_US_Mfg/ECN_2012_US_31SG1_with_ann.csv", clear names
    keep if naicsid == "31-33"
    keep company
    drop if _n == 2
    compress
    destring, replace
    gen iso3_o = "USA"
    gen country = "United States"
    gen variable = "Number of active employer enterprises"
    rename  company value
    gen year = `year'
    gen source_year = ${year}
    gen source = "2012 US Census"
    saveold "$ROOT/Data/Int/Census`year'_small", replace
    append using  "$ROOT/Data/Int/SDBS`year'"

    drop country
    rename value SDBS_N_ii
    label var SDBS_N_ii "Number of Firms"
    drop variable
    drop source_year
    rename source source_SDBS
    saveold "$ROOT/Data/Int/SDBS`year'_small", replace

  end

  * Uses the WB Data to find firm export propbilites by destination
  capture program drop make_WB_Enterprise
  program define make_WB_Enterprise
  syntax, year(integer)

    clear
    insheet using "$ROOT/Data/WB_Enterprise_Data/Name_crosswalk.tsv", clear names
    tempfile isoxw
    saveold `isoxw',replace

    clear
    insheet using "$ROOT/Data/WB_Enterprise_Data/Full_Trade.csv", clear names
    gen year_diff = abs(year - `year')
    gsort economy year_diff -year
    by economy: keep if _n == 1
    keep economy year percentoffirmsexportingdirectlya
    rename year source_year_ENT
    gen year = `year'
    rename economy country
    rename percentoffirmsexportingdirectlya export_probability
    label var export_probability "Export Probability from World Bank Enterprise Survey"
    merge 1:1 country using   `isoxw', keep(master match)
    assert _merge == 3
    drop _merge
    drop country
    saveold "$ROOT/Data/Int/WB_Enterprise_Data`year'_small", replace
  end



// Step 3
// Exporting Firm Data
capture program drop make_WB_EDD
program define make_WB_EDD
syntax, year(integer)

    * Load Total Exporter Statistics
    use "$ROOT/Data/EDD_Data/CY_manuf/CY_manuf.dta", clear
    gen EDD_N_iE = A1
    label var EDD_N_iE "Number of Total Exporters"
    // label var lN_iE "Log Number of Total Exporters"
    keep y c  EDD_N_iE
    rename y source_year_EDD
    rename c iso3_o
    saveold "$ROOT/Data/Int/EDD_CY_MAN", replace

    * Load Total Exporter Statistics
    use "$ROOT/Data/EDD_Data/CY_manuf/CY_manuf.dta", clear
    gen year_diff = abs(y - `year')
    gsort c year_diff -y
    by c: keep if _n == 1
    drop year_diff
    gen year = `year'
    rename y source_year_EDD
    rename c iso3_o

    tempfile EDD_CY_MAN
    save `EDD_CY_MAN', replace
    use "$ROOT/Data/EDD_Data/CYD_manuf/CYD_manuf.dta", clear

    rename y source_year_EDD
    rename c iso3_o
    merge m:1 iso3_o source_year_EDD using  `EDD_CY_MAN', keep(match)
    order year
    assert _merge == 3
    drop _merge
    gen source_EDD = "EDD"
    tab y
    gen EDD_X_ij = A1*A6i
    gen EDD_N_ij = (A1)
    rename d iso3_d
    keep year iso3* source_year_EDD A1 A6i A6ii source* EDD*
    saveold "$ROOT/Data/Int/WB_EDD`year'_small", replace

  end

  capture program drop make_OECD_TEC
  program define make_OECD_TEC
  syntax, year(integer)


    insheet using "$ROOT/Data/OECD/TES3_04052018220827756.csv", clear comma
    keep if sector == "C_D_E" // Keep manufacturers
    keep if itcs_flow ==  2 // For Exporter
    gen type =  cond(indicator==1,"N","X")
    keep reporter partner time value type
    drop if time == 2008
    tempfile temppre2008
    save `temppre2008', replace

    /* Import in Recent OECD Data */
      insheet using "$ROOT/Data/OECD/TEC3_REV4-en.csv", clear delim("|")
      keep if sector == "BCDE" // Keep manufacturers
      keep if flow ==  2 // For Exporter
      gen type =  cond(indicator==1,"N","X")
      drop  powercodecode referenceperiodcode flagcodes flow  unit indicator sector
      append using `temppre2008'

      preserve
        keep reporter time
        duplicates drop
        gen year_diff = abs(time - `year')
        gsort reporter year_diff -time
        by reporter: keep if _n == 1
        keep reporter time
        tempfile temp
        save `temp', replace
      restore

      merge m:1 reporter time using  `temp', keep(match)
      assert _merge == 3
      drop _merge

      // drop  powercodecode referenceperiodcode flagcodes flow  unit indicator sector
      reshape wide value, j(type) i(reporter partner time) string

      gen year = `year'
      rename time source_year_TEC

      rename valueN TEC_N_ij
      rename valueX TEC_X_ij
      rename  reporter  iso3_o
      rename  partner iso3_d
      saveold "$ROOT/Data/Int/OECDTEC_CDY1_MAN_`year'", replace

      keep if iso3_d == "TOTAL"
      drop iso3_d
      rename TEC_N_ij TEC_N_iE
      rename TEC_X_ij TEC_X_iE
      saveold "$ROOT/Data/Int/OECDTEC_CY1_MAN_`year'", replace

  end

  capture program drop make_OECD_TEC_ALL_FIRM
  program define make_OECD_TEC_ALL_FIRM
  syntax, year(integer)



      insheet using "$ROOT/Data/OECD/TES3_04052018220827756.csv", clear comma
      keep if sector == "C_D_E" // Keep manufacturers
      keep if itcs_flow ==  2 // For Exporter
      gen type =  cond(indicator==1,"N","X")
      keep reporter partner time value type
      drop if time == 2008
      tempfile temppre2008
      save `temppre2008', replace



    /* Import in Recent OECD Data */
    insheet using "$ROOT/Data/OECD/TEC3_REV4-en.csv", clear delim("|")
    keep if sector == "TOTAL" // Keep manufacturers
    keep if flow ==  2 // For Exporter
    gen type =  cond(indicator==1,"N","X")
    drop  powercodecode referenceperiodcode flagcodes flow  unit indicator sector
    append using `temppre2008'

    preserve
      keep reporter time
      duplicates drop
      gen year_diff = abs(time - `year')
      gsort reporter year_diff -time
      by reporter: keep if _n == 1
      keep reporter time
      tempfile temp
      save `temp', replace
    restore

    merge m:1 reporter time using  `temp', keep(match)
    assert _merge == 3
    drop _merge

    reshape wide value, j(type) i(reporter partner time) string

    gen year = `year'
    rename time source_year_TEC

    rename valueN TEC_TOTAL_N_ij
    rename valueX TEC_TOTAL_X_ij
    rename  reporter  iso3_o
    rename  partner iso3_d
    saveold "$ROOT/Data/Int/OECDTEC_CDY1_TOTAL_`year'", replace

    keep if iso3_d == "TOTAL"
    drop iso3_d
    rename TEC_TOTAL_N_ij TEC_TOTAL_N_iE
    rename TEC_TOTAL_X_ij TEC_TOTAL_X_iE
    saveold "$ROOT/Data/Int/OECDTEC_CY1_TOTAL_`year'", replace

  end


  capture program drop make_CHN_exporters
  program define make_CHN_exporters
  syntax, year(integer)





    // Gen N_ii
    import excel "$ROOT/Data/China_Statistical_Yearbook/Combined_Industrial_Enterprises_Indicators.xlsx", sheet("Sheet1") firstrow clear
    gen iso3_o = "CHN"
    rename NumberofEnterprisesunit CHN_N_ii
    rename Year year
    keep year CHN_N_ii iso3_o
    tempfile CHN1
    saveold `CHN1', replace

    // Gen X_ii
    import excel "$ROOT/Data/China_Statistical_Yearbook/Combined_Industrial_Enterprises_Indicators.xlsx", sheet("Sheet1") firstrow clear
    rename Year year
    keep if year == `year'
    gen iso3_o = "CHN"
    gen iso3_d = "CHN"
    rename NumberofEnterprisesunit CHN_N_ij
    rename BusinessRevenue CHN_X_ij_yuan
    replace CHN_X_ij_yuan = CHN_X_ij_yuan*10 //  Trying to line up with WIOD - FX rate
    keep year CHN_N_ij CHN_X_ij_yuan iso3_o iso3_d
    tempfile CHN2
    saveold `CHN2', replace


    // Put All china data together
    use "$ROOT/Data/China_Exports/Final/Final Data.dta", clear
    drop if isocode == "CHN"
    keep if year == `year'
    rename n CHN_N_ij
    rename value CHN_X_ij_USD
    replace CHN_X_ij_USD = CHN_X_ij_USD/1000000
    rename isocode iso3_d
    gen iso3_o = "CHN"
    cap collapse (sum) CHN_N_ij CHN_X_ij, by(iso3_o iso3_d year)
    append using `CHN2'
    merge m:1 iso3_o  year  using `CHN1' ,keep(master matched)


    assert _merge == 3
    drop _merge
    merge m:1 iso3_o  year  using "$ROOT/Data/Int/WIOD_FX" ,keep(master matched)
    assert _merge == 3
    drop _merge
    replace CHN_X_ij_USD = CHN_X_ij_yuan/US_FX if CHN_X_ij_USD == .
    drop US_FX
    saveold "$ROOT/DATA/Int/CHN_`year'", replace




  end

  capture program drop make_AUS_exporters
  program define make_AUS_exporters
  syntax, year(integer)

    insheet using "$ROOT/Data/Australia/AUS_cleaned.csv", clear names case
    drop ROW Dest_Country
    reshape long N_ij_ T_ij_ X_ij_, i(iso3_o  iso3_d) j(year)
    rename N_ij_ AUS_N_ij
    rename T_ij_ AUS_T_ij
    rename X_ij_ AUS_X_ij
    destring  AUS*, ignore(",-np") replace

    saveold "$ROOT/DATA/Int/AUS_all", replace

    gen year_diff = abs(year - `year')
    gsort iso3_d year_diff -year
    by iso3_d: keep if _n == 1
    merge m:1 iso3_o  year  using "$ROOT/Data/Int/WIOD_FX" ,keep(master matched)
    replace AUS_X_ij = AUS_X_ij/US_FX
    drop year year_diff _merge US_FX

    saveold "$ROOT/DATA/Int/AUS_`year'", replace

  end


  //Step 4
  //Gravity data
  capture program drop make_CEPII
  program define make_CEPII

    global CEPII_Dist "$ROOT/Data/CEPII/GeoDist/dist_cepii.dta"
    use $CEPII_Dist, clear
    rename iso_o iso3_o
    rename iso_d iso3_d
    replace iso3_o = "ROU" if iso3_o == "ROM"
    replace iso3_d = "ROU" if iso3_d == "ROM"
    saveold "$ROOT/Data/Int/CEPII_Dist_temp", replace

    global CEPII_Grav $ROOT/Data/CEPII/Gravity/gravdata_cepii.dta
    use $CEPII_Grav, clear
    replace iso3_o = "ROU" if iso3_o == "ROM"
    replace iso3_d = "ROU" if iso3_d == "ROM"
    saveold "$ROOT/Data/Int/CEPII_Grav_temp", replace

    global CEPII_Grav $ROOT/Data/CEPII/Gravity_dta_V202102/Gravity_V202102.dta
    use $CEPII_Grav, clear
    keep if year == 2012
    drop if distw == .
    drop if iso3_o == "ROW"
    drop if iso3_d == "ROW"
    drop if iso3_o == "TOT"

    // No Distance Data
    drop  if inlist(iso3_o, "ANT","CSK","DDR","SCG","SUN","YMD","YUG")
    drop  if inlist(iso3_d, "ANT","CSK","DDR","SCG","SUN","YMD","YUG")

    // No GDP Data
    drop if gdp_d == .
    drop if gdp_o == .

    encode iso3_o, gen(iso3_o_n)
    preserve
      keep iso3_o iso3_o_n
      rename iso3_o iso3_d
      rename iso3_o_n iso3_d_n
      duplicates drop
      tempfile x1
      saveold `x1', replace
    restore
    merge m:1 iso3_d using `x1', nogen

    gen iso3_o_n2 = iso3_o_n
    gen iso3_d_n2 = iso3_d_n
    order iso3_o_n2 iso3_d_n2     distw gdp_o gdp_d distw


    preserve
      keep iso3_o iso3_d iso3_o_n2 iso3_d_n2      distw gdp_o gdp_d distw tradeflow_baci pop_o pop_d
      order iso3_o iso3_d iso3_o_n2 iso3_d_n2     distw gdp_o gdp_d distw tradeflow_baci pop_o pop_d
      sort iso3_o iso3_d
      cap mkdir "$ROOT/Data/Int/DIST/"
      outsheet using "$ROOT/Data/Int/DIST/full.csv", replace  comma
    restore


    saveold "$ROOT/Data/Int/CEPII_Grav_full", replace

    keep iso3_o iso3_o_n2
    duplicates drop
    sort iso3_o_n2
    outsheet using "$ROOT/Data/Int/DIST/label.csv", replace  comma non




  end





//Step 5: survival probability data
capture program drop make_OECD_entry
program define make_OECD_entry
syntax, year(integer)
  global year = `year'

  insheet using "$ROOT/Data/OECD/SDBS Business Demography Indicators/SDBS_BDI_ISIC4_11042019204932050.csv", clear names case
  keep  if Variable == "2-year survival rate"
  keep if SEC == "10_33" // Keep MFG
  keep if SizeClass == "Total"
  keep if TIME == ${year}
  rename  LOCATION iso3_o
  rename  Value SurvivalRate2
  label var SurvivalRate2 "OECD 2-year survival rate"
  keep iso3_o SurvivalRate2
  drop if SurvivalRate2 >= 100
  saveold "$ROOT/DATA/Int/SurvivalRate2${year}", replace

  insheet using "$ROOT/Data/OECD/SDBS Business Demography Indicators/SDBS_BDI_ISIC4_11042019204932050.csv", clear names case
  keep  if Variable == "3-year survival rate"
  keep if SEC == "10_33" // Keep MFG
  keep if SizeClass == "Total"
  keep if TIME == ${year}
  rename  LOCATION iso3_o
  rename  Value SurvivalRate3
  label var SurvivalRate3 "OECD 2-year survival rate"
  keep iso3_o SurvivalRate3
  drop if SurvivalRate3 >= 100
  saveold "$ROOT/DATA/Int/SurvivalRate3${year}", replace

  insheet using "$ROOT/Data/OECD/SDBS Business Demography Indicators/SDBS_BDI_ISIC4_11042019204932050.csv", clear names case
  keep  if Variable == "1-year survival rate"
  keep if SEC == "10_33" // Keep MFG
  keep if SizeClass == "Total"
  keep if TIME == ${year}
  rename  LOCATION iso3_o
  rename  Value SurvivalRate1
  label var SurvivalRate1 "OECD 1-year survival rate"
  keep iso3_o SurvivalRate1
  drop if SurvivalRate1>= 100
  merge 1:1 iso3_o using "$ROOT/DATA/Int/SurvivalRate2${year}", nogen
  merge 1:1 iso3_o using "$ROOT/DATA/Int/SurvivalRate3${year}", nogen
  capture egen SurvivalRate1_average = mean(SurvivalRate1)
  capture label var SurvivalRate1_average "Average OECD 1-year survival rate"
  capture egen SurvivalRate2_average = mean(SurvivalRate2)
  capture label var SurvivalRate2_average "Average OECD 2-year survival rate"
  capture egen SurvivalRate3_average = mean(SurvivalRate3)
  capture label var SurvivalRate3_average "Average OECD 5-year survival rate"
  saveold "$ROOT/DATA/Int/SurvivalRate${year}", replace

end



// Step 6
// Stack and Construct Data
capture program drop stack_data
program define stack_data
syntax, year(integer)
  global year = `year'


  // Start with the Diagonals
  // global year = ${year}
  use "$ROOT/Data/Int/WIOT${year}_small", replace
  gen year = ${year}

  drop if inlist(iso3_o,"ROW","TOT")
  drop if inlist(iso3_d,"ROW","TOT")


  // Merge in BACI trade data
  merge 1:m iso3_o iso3_d year using "$ROOT/Data/Int/BACI${year}_small"
  drop _merge
  drop year
  

  // We Dont Use This
  /* 
  // Trade Cost Data
  gen year = ${year}
  merge 1:m iso3_o iso3_d year using "$ROOT/Data/Freight data/gravity_freight_data_unmerged", keep(master matched)
  keep if year == ${year}
  tab iso3_d _merge
  drop _merge
  drop year

  // Trade Cost Data
  cap gen year = ${year}
  merge 1:m iso3_o iso3_d year using "$ROOT/Data/Freight data/gravity_ESCAP_WB_TC_square" //, keep(master matched)
  keep if year == ${year}
  tab iso3_d _merge
  drop _merge
  drop year 
  */

  // Origin OECD Export Statistics
  merge m:1 iso3_o using "$ROOT/Data/Int/OECDTEC_CY1_MAN_${year}", keep(master matched)
  tab iso3_o _merge
  drop _merge

  // O-D OECD Survival Statistics
  merge m:1 iso3_o using "$ROOT/Data/Int/SurvivalRate${year}", keep(master matched)
  tab iso3_o _merge
  drop _merge

  // O-D OECD Export Statistics
  merge m:1 iso3_o iso3_d using "$ROOT/Data/Int/OECDTEC_CDY1_MAN_${year}", keep(master matched)
  tab iso3_o _merge
  drop _merge

  // Origin OECD Export Statistics
  merge m:1 iso3_o using "$ROOT/Data/Int/OECDTEC_CY1_TOTAL_${year}", keep(master matched)
  tab iso3_o _merge
  drop _merge

  // O-D OECD Export Statistics
  merge m:1 iso3_o iso3_d using "$ROOT/Data/Int/OECDTEC_CDY1_TOTAL_${year}", keep(master matched)
  tab iso3_o _merge
  drop _merge

  // O-D EDD Export Statistics
  merge 1:1 iso3_o iso3_d using "$ROOT/Data/Int/WB_EDD${year}_small", keep(master matched)
  tab iso3_o _merge
  drop _merge

  // Origin EDD cvExport Statistics
  merge m:1 iso3_o source_year_EDD using "$ROOT/Data/Int/EDD_CY_MAN", keep(master matched)
  tab iso3_o _merge
  drop _merge

  // Merge N_ii data OECD
  merge m:1 iso3_o  using "$ROOT/Data/Int/SDBS${year}_small", keep(master matched)
  tab iso3_o _merge
  drop _merge

  merge m:1 iso3_o  using "$ROOT/Data/Int/SSIS${year}_small", keep(master matched)
  tab iso3_o _merge
  drop _merge

  // Merge N_ii data OECD - WB
  merge m:1 iso3_o  using "$ROOT/Data/Int/WB_Enterprise_Data${year}_small", keep(master matched)
  tab iso3_o _merge
  drop _merge

  // AUS Data
  merge m:1 iso3_o iso3_d using "$ROOT/Data/Int/AUS_${year}", keep(master matched)
  tab iso3_o _merge
  drop _merge

  // CHN Data
  merge m:1 iso3_o iso3_d using "$ROOT/Data/Int/CHN_${year}", keep(master matched)
  tab iso3_o _merge
  drop _merge


  // Merge in Required Square Data
  merge m:1 iso3_o iso3_d using "$ROOT/Crosswalks/iso3_o_d"
  drop _merge
  cap gen year = ${year}


  merge m:1 iso3_o iso3_d using "$ROOT/Data/Int/CEPII_Dist_temp", keep(master matched)
  tab iso3_o _merge
  drop _merge

  replace year = ${year}
  merge m:1 iso3_o iso3_d year using "$ROOT/Data/Int/CEPII_Grav_temp", keep(master matched)
  tab iso3_o _merge
  drop _merge
  count


  // Distribute imputed survival rates
  // To Do: get better survival rate data
  replace SurvivalRate1_average = 85 if SurvivalRate1_average == . // For years without data
  replace SurvivalRate2_average = 75  if SurvivalRate2_average == . // For years without data
  replace SurvivalRate3_average = 65   if SurvivalRate3_average == . // For years without data

  forv i=1(1)3 {
    di "Suvival rate `i'"
    capture gen SurvivalRate`i'= .
    egen SurvivalRate`i'_average_ = mean(SurvivalRate`i'_average)
    replace   SurvivalRate`i'_average = SurvivalRate`i'_average_
    drop SurvivalRate`i'_average_
    gen SurvivalRate`i'_combined = SurvivalRate`i'
    replace SurvivalRate`i'_combined = SurvivalRate`i'_average if SurvivalRate`i'_combined ==.
    replace SurvivalRate`i'_combined = SurvivalRate`i'_combined/100

  }

  // Use BACI for diag (as a starting point)
  gen X_ij = BACI_v if iso3_o!= iso3_d
  gen X_ij_source = "BACI" if X_ij != .

  replace X_ij=WIOD_X_ij if X_ij == .
  replace X_ij_source="WIOD" if !missing(X_ij) & X_ij_source == ""

  // Create own firm data
  cap drop N_ii*
  gen N_ii = SDBS_N_ii
  gen N_ii_source = "SDBS" if !missing(N_ii)

  replace N_ii = SSIS_N_ii if missing(N_ii) // ALternative Year Called SSIS
  replace N_ii_source = "SDBS" if !missing(N_ii) & N_ii_source == ""

  gen N_ii_EDD_impute = EDD_N_iE/(export_probability/100)
  replace N_ii = N_ii_EDD_impute if missing(N_ii)
  replace N_ii_source = "EDD/WBES" if !missing(N_ii) & N_ii_source == ""

  replace N_ii = AUS_N_ij if iso3_o == iso3_d & missing(N_ii)
  replace N_ii_source = "AUS" if !missing(N_ii) & N_ii_source == ""

  replace N_ii = CHN_N_ii if iso3_o == "CHN"
  replace N_ii_source = "CHN" if !missing(N_ii) & N_ii_source == ""

  sort iso3_o -N_ii
  by iso3_o: replace N_ii = N_ii[1] if missing(N_ii)
  by iso3_o: replace N_ii_source = N_ii_source[1] if !missing(N_ii) & N_ii_source == ""

  tab N_ii_source

  // Create exporter data
  gen N_ij = TEC_N_ij
  gen N_ij_source = "TEC" if !missing(N_ij)

  replace N_ij = TEC_TOTAL_N_ij if missing(N_ij) & iso3_o == "KOR"
  replace N_ij_source = "TEC*" if !missing(N_ij) & N_ij_source == ""

  replace N_ij = EDD_N_ij if missing(N_ij)
  replace N_ij_source = "EDD" if !missing(N_ij) & N_ij_source == ""

  replace N_ij = AUS_N_ij if missing(N_ij)
  replace N_ij_source = "AUS" if !missing(N_ij) & N_ij_source == ""

  replace N_ij = CHN_N_ij if missing(N_ij)
  replace N_ij_source = "CHN" if !missing(N_ij) & N_ij_source == ""

  replace N_ij = N_ii if iso3_o == iso3_d
  replace N_ij_source = N_ii_source if !missing(N_ij) & N_ij_source == ""

  tab N_ij_source

  // Create x_bar data from consistent source
  cap drop x_bar
  gen x_bar = TEC_X_ij/TEC_N_ij
  gen x_bar_source = "TEC" if !missing(x_bar) 

  order x_bar
  sort x_bar
  replace x_bar = TEC_TOTAL_X_ij/TEC_TOTAL_N_ij if missing(x_bar) & iso3_o == "KOR"
  replace x_bar_source = "TEC" if !missing(x_bar) & x_bar_source == ""

  replace x_bar = A6i/1000000 if missing(x_bar)
  replace x_bar_source = "EDD" if !missing(x_bar) & x_bar_source == ""

  replace x_bar = AUS_X_ij/AUS_N_ij if missing(x_bar)
  replace x_bar_source = "AUS" if !missing(x_bar) & x_bar_source == ""

  replace x_bar = X_ij/AUS_N_ij if missing(x_bar)
  replace x_bar_source = "AUS" if !missing(x_bar) & x_bar_source == ""

  replace x_bar = CHN_X_ij_USD/CHN_N_ij if iso3_o == "CHN"
  replace x_bar_source = "CHN" if !missing(x_bar) & x_bar_source == ""

  // There is Censoring in the EDD data
  replace x_bar = X_ij/N_ij if iso3_o == iso3_d  & missing(x_bar)
  replace x_bar_source = X_ij_source +" + "+N_ij_source if !missing(x_bar) & x_bar_source == ""

  // DO we need to adjust the X_ii for the diagonals
  // Rodrigo, 10:20 AM
  // To merge the datasets, I would jsut scale X_ii such that the ratio of total exports in WIOD and our data
  // So that the trade matrix is all cosnsitent
  // for each country

  // Use a Consistent X_ij
  replace X_ij = x_bar*N_ij
  replace X_ij = BACI_v if X_ij == .
  replace X_ij = WIOD_X_ij if X_ij == .
  replace X_ij_source = x_bar_source if !missing(X_ij) &  X_ij_source == ""

  replace x_bar = X_ij/N_ij if missing(x_bar)
  replace x_bar_source = X_ij_source +" + "+N_ij_source if !missing(x_bar) & x_bar_source == ""

  // Diagnostics
  // gen X_ij_1 = x_bar*N_ij
  // gen lX_ij_1 = log(X_ij_1)
  // gen lWIOD_X_ij = log(WIOD_X_ij)
  // gen lBACI_v = log(BACI_v)
  // corr lX_ij_1 lWIOD_X_ij lBACI_v
  // regress lX_ij_1 lWIOD_X_ij
  // order lX_ij_1 lWIOD_X_ij X_ij_1 WIOD_X_ij
  // br if iso3_o == "CHN" || iso3_o == "DEU"

  gen n_ij = N_ij/N_ii*SurvivalRate1_combined
  gen n_ij_survival2 = N_ij/N_ii*SurvivalRate2_combined
  gen n_ij_survival3 = N_ij/N_ii*SurvivalRate3_combined
  gen n_ij_allsurvival = N_ij/N_ii

  order iso* X_ij N_ii N_ij x_bar n_ij
  sort iso3_o iso3_d

  gen ldist = log(dist)
  gen ldistw = log(distw)
  gen ln_ij =log(n_ij)
  gen lx_bar=log(x_bar)
  gen lN_ii = log(N_ii)
  gen lN_ij = log(N_ij)
  encode iso3_o, gen(o)
  encode iso3_d, gen(d)

  count if missing(N_ij)
  tab iso3_o if missing(N_ij)
  count

  saveold "$ROOT/Data/Int/stack_data_${year}_B", replace

end




capture program drop STACK_ALL_YEARS
program define STACK_ALL_YEARS

  clear
  forvalues y = 2000(1)2014 {
    di "`y'"
    global year = `y'
    append using "$ROOT/Data/Int/stack_data_${year}_B"
  }

  drop if inlist(iso3_o,"ROW","TOT")
  gen lX_ij = log(X_ij)
  egen group_oy = group(iso3_o year)
  egen group_dy = group(iso3_d year)

  gen gdpcap_d_2000_ = gdpcap_d if year == 2000
  gen gdpcap_o_2000_ = gdpcap_o if year == 2000
  bys iso3_d: egen gdpcap_d_2000 = max(gdpcap_d_2000_)
  bys iso3_o: egen gdpcap_o_2000 = max(gdpcap_o_2000_)
  drop gdpcap_d_2000_ gdpcap_o_2000_

  // Use $9,000 as cutoff 
  gen rich_dest = cond(gdpcap_d_2000>9000,1,0)
  gen rich_orig = cond(gdpcap_o_2000>9000,1,0)

  // preserve
  //   keep if year == 2000
  //   keep iso3_d gdpcap_d
  //   duplicates drop
  //   summ gdpcap_d, d
  //   di "`r(p50)'"
  //   global CUTOFF  `r(p50)'
  // restore
  // gen rich_dest = cond(gdpcap_d_2000>$CUTOFF,1,0)
  // gen rich_orig = cond(gdpcap_o_2000>$CUTOFF,1,0)
  // bys iso3_o iso3_d: egen count = count(lx_bar)
  // tab count
  // keep if count == 15

  label var lx_bar "ln(x_ij_bar)"
  label var lN_ij "ln(N_ij)"
  label var lX_ij "ln(X_ij)"
  label var ldistw "ln(Distw_ij)"
  saveold "$ROOT/Data/Int/stack_data_B_", replace

  use "$ROOT/Data/Int/stack_data_B_", clear
  merge 1:1 iso3_o iso3_d year using "$ROOT/Data/TRAINS/collapsed/T_all_collapse.dta", keep(master matched)
  replace simpleAHS_w = 0 if iso3_o == iso3_d
  replace simpleAHS_uw = 0 if iso3_o == iso3_d
  drop _merge

  merge 1:1 iso3_o iso3_d year using "$ROOT/Data/TRAINS/collapsed/Teti_all_collapse.dta", keep(master matched)
  replace TsimpleAHS_uw = 0 if iso3_o == iso3_d
  drop _merge

  saveold "$ROOT/Data/Int/stack_data_B", replace

end


capture program drop GSP_2012_data
program define GSP_2012_data

  import excel "$ROOT/Data/GSP/gspmfnrates/mfn_gsp_2012.xlsx", sheet("gsp_dummy_2012") cellrange(A1:J200) firstrow clear
  rename exporter_code iso3_o
  rename (*_gsp) (gsp_ind*)
  saveold "$ROOT/Data/Int/gsp_dummy_2012", replace

  import excel "$ROOT/Data/GSP/gspmfnrates/mfn_gsp_2012.xlsx", sheet("gsp_rates_2012") cellrange(A1:J200) firstrow clear
  rename exporter_code iso3_o
  rename (*_gsp_rate) (gsp_rate*)
  saveold "$ROOT/Data/Int/gsp_rates_2012", replace

  import excel "$ROOT/Data/GSP/gspmfnrates/mfn_gsp_2012.xlsx", sheet("mfn_rates_2012") cellrange(A1:J200) firstrow clear
  rename exporter_code iso3_o
  rename (*_mfn_rate) (mfn_rate*)
  saveold "$ROOT/Data/Int/mfn_rates_2012", replace

  import excel "$ROOT/Data/GSP/gspmfnrates/mfn_gsp_2012.xlsx", sheet("gsp_rates_2012_wt") cellrange(A1:J200) firstrow clear
  rename exporter_code iso3_o
  rename (*_gsp_rate) (gsp_rate_wt*)

  saveold "$ROOT/Data/Int/gsp_rates_2012_wt", replace

  import excel "$ROOT/Data/GSP/gspmfnrates/mfn_gsp_2012.xlsx", sheet("mfn_rates_2012_wt") cellrange(A1:J200) firstrow clear
  rename exporter_code iso3_o
  rename (*_mfn_rate) (mfn_rate_wt*)
  saveold "$ROOT/Data/Int/mfn_rates_2012_wt", replace


  clear
  use "$ROOT/Data/Int/gsp_dummy_2012"
  merge 1:1 iso3_o using "$ROOT/Data/Int/gsp_rates_2012", nogen
  merge 1:1 iso3_o using "$ROOT/Data/Int/mfn_rates_2012", nogen
  merge 1:1 iso3_o using "$ROOT/Data/Int/gsp_rates_2012_wt", nogen
  merge 1:1 iso3_o using "$ROOT/Data/Int/mfn_rates_2012_wt", nogen



  reshape long gsp_ind gsp_rate mfn_rate gsp_rate_wt mfn_rate_wt, i(iso3_o) j(iso3_d) string
  destring  gsp_ind gsp_rate* mfn_rate*, force replace

  // Generate Ratios
  // Cap at 1
  gen delta = gsp_ind*(100+mfn_rate)/(100+gsp_rate)
  replace delta = 1 if gsp_ind == 0
  replace delta = 1 if delta < 1

  gen delta_ind = gsp_ind*(100+mfn_rate)/100
  replace delta_ind = 1 if delta_ind == 0
  replace delta_ind = 1 if delta_ind < 1

  gen delta_wt = gsp_ind*(100+mfn_rate_wt)/(100+gsp_rate_wt)
  replace delta_wt = 1 if gsp_ind == 0
  replace delta_wt = 1 if delta_wt < 1

  drop if iso3_o == iso3_d
  rename iso3_d iso3_d2 
  saveold "$ROOT/Data/Int/GSP_2012_data", replace

end


capture program drop DATA_EXPORT_AGG
program define DATA_EXPORT_AGG
    version 13
    syntax [namelist], year(integer)

  global year = `year'
   // global year = 2012
  // use "$ROOT/Data/Int/stack_data_${year}", clear


  use  "$ROOT/Data/Int/stack_data_B", clear
  drop o d
  egen group = group(iso3_o iso3_d)
  tsset group year
  gen simpleAHS_w_l5 = l5.simpleAHS_w 
  gen simpleMFN_w_l5 = l5.simpleMFN_w
  gen simpleAHS_uw_l5 = l5.simpleAHS_uw 
  gen simpleMFN_uw_l5 = l5.simpleMFN_uw
  gen simpleAHS_uw_l10 = l10.simpleAHS_uw 
  gen simpleMFN_uw_l10 = l10.simpleMFN_uw


  gen fta_wto_l5     = l5.fta_wto
  gen gatt_o_l5      = l5.gatt_o
  gen gatt_d_l5      = l5.gatt_d
  gen fta_wto_l10     = l10.fta_wto
  gen gatt_o_l10      = l10.gatt_o
  gen gatt_d_l10      = l10.gatt_d

  keep if year == ${year}

  // Summary Table: Source Simple
  preserve
      drop if iso3_o == iso3_d
      gcollapse (sum) BACI_v // X_ij WIOD_X_ij
      gen sample = "Aggregate"
      saveold  "$ROOT/Data/Int/Aggregate_Flows", replace
  restore


  preserve
      drop if iso3_o == iso3_d
      gen tariff_use = simpleAHS_w/100
      gen ltariff_use = log(1+tariff_use)
      gen ltariff_use2 = ltariff_use
      gen ldist_use = ldistw
      order ln_ij lx_bar ldistw ltariff_use

      gen valid = 1
      replace valid = 0 if ln_ij == .
      replace valid = 0 if lx_bar == .
      replace valid = 0 if ldistw == .
      replace valid = 0 if ltariff_use == .
      replace valid = 0 if lx_bar <= -20
      gen flows = valid*X_ij
      count

      gcollapse (sum) flows  X_ij WIOD_X_ij valid, by(iso3_o  year) 
      gsort -X_ij -valid
      gen include = ""
      drop if X_ij == 0
      replace include = "to include" if flows > 0
      replace include = "to include" if WIOD_X_ij != 0
      gsort -include -X_ij -valid
      drop flows
      outsheet using "$ROOT/Data/Int/WIOD_sampleB/${year}_CF_sample.csv", replace comma
  restore

  drop if distw == .
  preserve
    keep if iso3_o== iso3_d
    keep iso3_o iso3_d  X_ij N_ii N_ij   *source SurvivalRate1*
    rename SurvivalRate1 SurvivalRate_OECD
    rename SurvivalRate1_average SurvivalRate_imputed
    drop SurvivalRate1_combined
    gen N = _n
    saveold  "$ROOT/Data/Int/temp", replace
    collapse (count) count = SurvivalRate_OECD count_imputed = SurvivalRate_imputed
    gen label = "Countries with Survival Rate Data"
    list
    outsheet using "$ROOT/Data/Int/WIOD_sampleB/${year}_Survival_Rates.csv", replace comma
  restore

  // Drop crazy outliers
  drop if lx_bar <= -20
  reghdfe ln_ij lx_bar ldistw, absorb(iso3_o iso3_d)
  count
  drop if !e(sample)
  count

  // Source Countries
  preserve
    keep iso3_o
    duplicates drop
    count
  restore

  // Destination Countries
  preserve
    keep iso3_d
    duplicates drop
    count
  restore

  // Fill In blanks for own.
  replace contig = 1 if iso3_o == iso3_d
  tab contig, missing

  replace fta_wto = 0 if fta_wto == .
  replace fta_wto = 1 if iso3_o == iso3_d
  tab fta_wto, missing

  gen dist_own  = distw if  iso3_o == iso3_d
  bys iso3_o: egen dist_ii = min(dist_own)
  bys iso3_d: egen dist_jj = min(dist_own)
  drop dist_own

  // Dummy Variables (depreciated)
  gen dummy1 = 0
  gen dummy2 = 0
  gen dummy3 = 0
  gen dummy4 = 0



  // Unweighted Tariffs
  preserve
    gen tariff_use = simpleAHS_uw/100
    gen ltariff_use = log(1+tariff_use)
    gen ltariff_use2 = ltariff_use
    gen ldist_use = ldistw

    reghdfe ln_ij lx_bar ldistw ltariff_use, absorb(iso3_o iso3_d)
    count
    drop if !e(sample)
    count
    reghdfe lX_ij  ltariff_use ldistw, absorb(iso3_o iso3_d)
    count
    drop if !e(sample)
    count
    

    encode iso3_o, gen(oN)
    encode iso3_d, gen(dN)
    summ iso3_o iso3_d

    assert  ldistw != .
    keep  oN dN lN_ij lx_bar ltariff_use gdpcap_o_2 gdpcap_d_2  rich_orig rich_dest ln_ij gdp_o gdp_d colony contig fta_wto comcur comlang_ethno  dummy1 dummy2 dummy3 dummy4 ldist_use ltariff_use2 
    order oN dN lN_ij lx_bar ltariff_use gdpcap_o_2 gdpcap_d_2  rich_orig rich_dest ln_ij gdp_o gdp_d colony contig fta_wto comcur comlang_ethno  dummy1 dummy2 dummy3 dummy4 ldist_use ltariff_use2 
    cap mkdir "$ROOT/Data/Int/WIOD_sampleB/"
    outsheet using "$ROOT/Data/Int/WIOD_sampleB/${year}_Tuw.csv", replace non nol comma
    outsheet using "$ROOT/Data/Int/WIOD_sampleB/${year}_Tuw_lab.csv", replace 
  restore


  // Unweighted Tariffs - Teti
  preserve
    gen tariff_use = TsimpleAHS_uw/100
    gen ltariff_use = log(1+tariff_use)
    gen ltariff_use2 = ltariff_use
    gen ldist_use = ldistw

    reghdfe ln_ij lx_bar ldistw ltariff_use, absorb(iso3_o iso3_d)
    count
    drop if !e(sample)
    count
    reghdfe lX_ij  ltariff_use ldistw, absorb(iso3_o iso3_d)
    count
    drop if !e(sample)
    count
    

    encode iso3_o, gen(oN)
    encode iso3_d, gen(dN)

    assert  ldistw != .
    keep  oN dN lN_ij lx_bar ltariff_use gdpcap_o_2 gdpcap_d_2  rich_orig rich_dest ln_ij gdp_o gdp_d colony contig fta_wto comcur comlang_ethno  dummy1 dummy2 dummy3 dummy4 ldist_use ltariff_use2 
    order oN dN lN_ij lx_bar ltariff_use gdpcap_o_2 gdpcap_d_2  rich_orig rich_dest ln_ij gdp_o gdp_d colony contig fta_wto comcur comlang_ethno  dummy1 dummy2 dummy3 dummy4 ldist_use ltariff_use2 
    cap mkdir "$ROOT/Data/Int/WIOD_sampleB/"
    outsheet using "$ROOT/Data/Int/WIOD_sampleB/${year}_Tetiuw.csv", replace non nol comma
    outsheet using "$ROOT/Data/Int/WIOD_sampleB/${year}_Tetiuw_lab.csv", replace 
  restore


















  if $year != 2012 {
    exit
  }
  else {
  }

  // Unweighted Tariffs - Teti - China Origin
  preserve
    replace rich_orig = cond(iso3_o=="CHN",1,0)
    replace rich_dest = cond(iso3_d=="CHN",1,0)

    gen tariff_use = TsimpleAHS_uw/100
    gen ltariff_use = log(1+tariff_use)
    gen ltariff_use2 = ltariff_use
    gen ldist_use = ldistw

    reghdfe ln_ij lx_bar ldistw ltariff_use, absorb(iso3_o iso3_d rich_orig rich_dest)
    count
    drop if !e(sample)
    count
    reghdfe lX_ij  ltariff_use ldistw, absorb(iso3_o iso3_d rich_orig rich_dest)
    count
    drop if !e(sample)
    count
    

    encode iso3_o, gen(oN)
    encode iso3_d, gen(dN)

    assert  ldistw != .
    keep  oN dN lN_ij lx_bar ltariff_use gdpcap_o_2 gdpcap_d_2  rich_orig rich_dest ln_ij gdp_o gdp_d colony contig fta_wto comcur comlang_ethno  dummy1 dummy2 dummy3 dummy4 ldist_use ltariff_use2 
    order oN dN lN_ij lx_bar ltariff_use gdpcap_o_2 gdpcap_d_2  rich_orig rich_dest ln_ij gdp_o gdp_d colony contig fta_wto comcur comlang_ethno  dummy1 dummy2 dummy3 dummy4 ldist_use ltariff_use2 
    outsheet using "$ROOT/Data/Int/WIOD_sampleB/${year}_Tetiuw_China.csv", replace non nol comma
  restore




  // Sample
  preserve
    gen tariff_use = simpleAHS_uw/100
    gen ltariff_use = log(1+tariff_use)
    gen ltariff_use2 = ltariff_use
    gen ldist_use = ldistw

    reghdfe ln_ij lx_bar ldistw ltariff_use, absorb(iso3_o iso3_d)
    count
    drop if !e(sample)
    count
    reghdfe lX_ij  ltariff_use ldistw, absorb(iso3_o iso3_d)
    count
    drop if !e(sample)
    count
    
    keep iso3_o
    duplicates drop
    save  "$ROOT/Data/Int/estimation_sample_${year}", replace
  restore


  // Summary Table: Source Simple
  preserve
      drop if iso3_o == iso3_d

      gen ltariff_use = log(1+simpleAHS_uw/100)

      reghdfe ln_ij lx_bar ldistw ltariff_use, absorb(iso3_o iso3_d)
      count
      drop if !e(sample)
      count
      reghdfe lX_ij  ltariff_use ldistw, absorb(iso3_o iso3_d)
      count
      drop if !e(sample)
      count
      
      gcollapse (sum) BACI_v  
      gen sample = "In Sample"
      append using "$ROOT/Data/Int/Aggregate_Flows"
      list 
      outsheet using "$ROOT/Output/${year}_Aggregate_Flows.csv", replace 
  restore





  // Summary Table: Source Full
  preserve
    gen ltariff_use = log(1+simpleAHS_uw/100)
    reghdfe ln_ij lx_bar ldistw ltariff_use, absorb(iso3_o iso3_d)
    count
    drop if !e(sample)
    count
    drop X*source
    collapse (count) count = ln_ij , by(iso3_o  rich_orig *source  )
    gsort iso3_o -count
    outsheet using "$ROOT/Output/${year}_Full_sources.csv", replace comma
    texsave  using "$ROOT/Output/${year}_Full_sources.tex", replace 
  restore

  // Summary Table: Statistics
  preserve
    gen ltariff_use = log(1+simpleAHS_uw/100)
    reghdfe ln_ij lx_bar ldistw ltariff_use, absorb(iso3_o iso3_d)
    count
    drop if !e(sample)
    count
    drop X*source
    collapse (count) count = ln_ij  (mean) ln_ij_mean = ln_ij lx_bar_mean = lx_bar (sd) ln_ij_sd = ln_ij lx_bar_sd = lx_bar , by(iso3_o  rich_orig  )
    gsort iso3_o -count
    format ln_ij_mean lx_bar_mean ln_ij_sd lx_bar_sd %12.2fc
    replace ln_ij_mean = round(ln_ij_mean,.01)
    replace lx_bar_mean = round(lx_bar_mean,.01)
    replace ln_ij_sd = round(ln_ij_sd,.01)
    replace lx_bar_sd = round(lx_bar_sd,.01)
    outsheet using "$ROOT/Output/${year}_Sample_summary.csv", replace comma
    texsave  using "$ROOT/Output/${year}_Sample_summary.tex", replace 
  restore

  // Main Sample - Weighted Tariffs
  preserve
    gen tariff_use = simpleAHS_w/100
    gen ltariff_use = log(1+tariff_use)
    gen ltariff_use2 = ltariff_use
    gen ldist_use = ldistw

    reghdfe ln_ij lx_bar ldistw ltariff_use, absorb(iso3_o iso3_d)
    drop if !e(sample)
    count

    encode iso3_o, gen(oN)
    encode iso3_d, gen(dN)

    assert  ldistw != .
    keep  oN dN lN_ij lx_bar ltariff_use gdpcap_o_2 gdpcap_d_2  rich_orig rich_dest ln_ij gdp_o gdp_d colony contig fta_wto comcur comlang_ethno  dummy1 dummy2 dummy3 dummy4 ldist_use ltariff_use2 
    order oN dN lN_ij lx_bar ltariff_use gdpcap_o_2 gdpcap_d_2  rich_orig rich_dest ln_ij gdp_o gdp_d colony contig fta_wto comcur comlang_ethno  dummy1 dummy2 dummy3 dummy4 ldist_use ltariff_use2 
    cap mkdir $ROOT/Data/Int/WIOD_sampleB/
    outsheet using "$ROOT/Data/Int/WIOD_sampleB/${year}_T.csv", replace non nol comma
    /* outsheet using "$ROOT/Data/Int/WIOD_sampleB/${year}_T_lab.csv", replace  */
  restore


  // Unweighted Tariffs IV
  preserve
    gen ldist_use = ldistw
    gen ltariff_IV = log(1+(simpleMFN_uw)/100)*(1-fta_wto)*gatt_o*gatt_d
    gen ltariff_IV2 = ltariff_IV 

    reghdfe ln_ij lx_bar ldistw ltariff_IV, absorb(iso3_o iso3_d)
    drop if !e(sample)

    encode iso3_o, gen(oN)
    encode iso3_d, gen(dN)

    assert  ldistw != .
    keep  oN dN lN_ij lx_bar ltariff_IV gdpcap_o_2 gdpcap_d_2  rich_orig rich_dest ln_ij gdp_o gdp_d colony contig fta_wto comcur comlang_ethno  dummy1 dummy2 dummy3 dummy4 ldist_use ltariff_IV2
    order oN dN lN_ij lx_bar ltariff_IV gdpcap_o_2 gdpcap_d_2  rich_orig rich_dest ln_ij gdp_o gdp_d colony contig fta_wto comcur comlang_ethno  dummy1 dummy2 dummy3 dummy4 ldist_use ltariff_IV2 
    cap mkdir $ROOT/Data/Int/WIOD_sampleB/
    outsheet using $ROOT/Data/Int/WIOD_sampleB/${year}_TuwIV.csv, replace non nol comma
  restore


  
  // Unweighted Tariffs IV - first difference - 10 year
  preserve
    gen ldist_use = ldistw
    gen ltariff_IV = (log(1+(simpleMFN_uw)/100) - log(1+(simpleMFN_uw_l10)/100))*(1-fta_wto_l10)*gatt_o_l10*gatt_d_l10
    gen ltariff_IV2 = ltariff_IV 


    reghdfe ln_ij lx_bar ldistw ltariff_IV, absorb(iso3_o iso3_d)
    drop if !e(sample)

    reghdfe ln_ij lx_bar ldistw ltariff_IV, absorb(iso3_o iso3_d)
    drop if !e(sample)

    reghdfe lX_ij  ltariff_IV ldistw, absorb(iso3_o iso3_d)
    count
    drop if !e(sample)

    count
    
    encode iso3_o, gen(oN)
    encode iso3_d, gen(dN)
    assert  ldistw != .
    keep  oN dN lN_ij lx_bar ltariff_IV gdpcap_o_2 gdpcap_d_2  rich_orig rich_dest ln_ij gdp_o gdp_d colony contig fta_wto comcur comlang_ethno  dummy1 dummy2 dummy3 dummy4 ldist_use ltariff_IV2
    order oN dN lN_ij lx_bar ltariff_IV gdpcap_o_2 gdpcap_d_2  rich_orig rich_dest ln_ij gdp_o gdp_d colony contig fta_wto comcur comlang_ethno  dummy1 dummy2 dummy3 dummy4 ldist_use ltariff_IV2 
    cap mkdir "$ROOT/Data/Int/WIOD_sampleB/"
    outsheet using "$ROOT/Data/Int/WIOD_sampleB/${year}_TuwIV2.csv", replace non nol comma
  restore
  

  // Data Histogram
  preserve
    label var ln_ij "Log of Exporter Firm Share"
    hist ln_ij, xlabel(-9.2 "0.01%" -6.9 "0.1%"    -4.6 	"1%"    -2.30 	"10%"   0 	"100%")
    graph export "$ROOT/Output/Hist_${year}B.pdf", replace
  restore

  gen tariff_use = simpleAHS_uw/100
  gen ltariff_use = log(1+tariff_use)
  gen ltariff_use2 = ltariff_use
  gen ldist_use = ldistw
    
  //Export estimation sample without own trade info
  preserve
    
    drop if iso3_o== iso3_d
    encode iso3_o, gen(oN)
    encode iso3_d, gen(dN)
    keep  oN dN lN_ij lx_bar ldistw gdpcap_o_2 gdpcap_d_2  rich_orig rich_dest ln_ij gdp_o gdp_d colony contig fta_wto  comcur comlang_ethno  dummy1 dummy2 dummy3 dummy4 ldist_use ltariff_use2 
    order oN dN lN_ij lx_bar ldistw gdpcap_o_2 gdpcap_d_2  rich_orig rich_dest ln_ij gdp_o gdp_d colony contig fta_wto  comcur comlang_ethno  dummy1 dummy2 dummy3 dummy4 ldist_use ltariff_use2 
    outsheet using "$ROOT/Data/Int/WIOD_sampleB/${year}_nonii_.csv", replace non nol comma
  restore

  //Export estimation sample without imputed survival rate
  preserve
    
    drop if SurvivalRate1== .
    encode iso3_o, gen(oN)
    encode iso3_d, gen(dN)
    keep  oN dN lN_ij lx_bar ldistw gdpcap_o_2 gdpcap_d_2  rich_orig rich_dest ln_ij gdp_o gdp_d colony contig fta_wto  comcur comlang_ethno  dummy1 dummy2 dummy3 dummy4 ldist_use ltariff_use2 
    order oN dN lN_ij lx_bar ldistw gdpcap_o_2 gdpcap_d_2  rich_orig rich_dest ln_ij gdp_o gdp_d colony contig fta_wto  comcur comlang_ethno  dummy1 dummy2 dummy3 dummy4 ldist_use ltariff_use2 
    outsheet using "$ROOT/Data/Int/WIOD_sampleB/${year}_survival_.csv", replace non nol comma
  restore

  //Export estimation sample with 3-year survival rate
  preserve
    
    gen ln_ij_survival3 = log(n_ij_survival3)
    encode iso3_o, gen(oN)
    encode iso3_d, gen(dN)
    keep  oN dN lN_ij lx_bar ldistw gdpcap_o_2 gdpcap_d_2  rich_orig rich_dest ln_ij_survival3 gdp_o gdp_d colony contig fta_wto comcur comlang_ethno  dummy1 dummy2 dummy3 dummy4 ldist_use ltariff_use2 
    order oN dN lN_ij lx_bar ldistw gdpcap_o_2 gdpcap_d_2  rich_orig rich_dest ln_ij_survival3 gdp_o gdp_d colony contig fta_wto comcur comlang_ethno  dummy1 dummy2 dummy3 dummy4 ldist_use ltariff_use2 
    outsheet using "$ROOT/Data/Int/WIOD_sampleB/${year}_survival_alt3_.csv", replace non nol comma
  restore

  //Export estimation sample under the assumption of n_ii = 1
  preserve
    
    gen ln_ij_allsurvival = log(n_ij_allsurvival)
    encode iso3_o, gen(oN)
    encode iso3_d, gen(dN)
    keep  oN dN lN_ij lx_bar ldistw gdpcap_o_2 gdpcap_d_2  rich_orig rich_dest ln_ij_allsurvival gdp_o gdp_d colony contig fta_wto comcur comlang_ethno  dummy1 dummy2 dummy3 dummy4 ldist_use ltariff_use2 
    order oN dN lN_ij lx_bar ldistw gdpcap_o_2 gdpcap_d_2  rich_orig rich_dest ln_ij_allsurvival gdp_o gdp_d colony contig fta_wto comcur comlang_ethno  dummy1 dummy2 dummy3 dummy4 ldist_use ltariff_use2 
    outsheet using "$ROOT/Data/Int/WIOD_sampleB/${year}_survival_all_.csv", replace non nol comma
  restore

end



*Step 7: Export data for counterfactual analysis
capture program drop DATA_EXPORT_CF_BIG
program define DATA_EXPORT_CF_BIG
  version 13
  syntax [namelist], year(integer) 


    global year = `year'



    // Get GSP Data
    use "$ROOT/Data/Int/stack_data_B", clear
    keep if year == 2012

    bys iso3_d: egen EU_D = total(eu_d) 
    replace eu_d = 1 if EU_D>0
    replace eu_d = 0 if eu_d == .

    bys iso3_o: egen EU_O = total(eu_o) 
    replace eu_o = 1 if EU_O>0 
    replace eu_o = 0 if eu_o == .

    // Make Sure Rich -> Rich isn't GSP
    replace gsp_d_d = 0 if rich_dest == rich_orig
    replace gsp_d_d = 0 if rich_orig == 1
    replace gsp_d_d = 0 if eu_o == 1
    replace gsp_d_d = 0 if  eu_o == 1 & eu_d == 1


    replace acp_to_eu = 0 if rich_dest == rich_orig
    replace acp_to_eu = 0 if rich_orig == 1
    replace acp_to_eu = 0 if eu_o == 1
    replace acp_to_eu = 0 if  eu_o == 1 & eu_d == 1

    tab   iso3_o gsp_d_d

    // Make Sure GSP from EU countries is unified
    replace acp_to_eu = 1 if eu_d == 1 & gsp_d_d == 1
    bys iso3_o: egen EU_GSP = total(acp_to_eu) 
    replace gsp_d_d = 1 if EU_GSP > 5 & eu_d == 1 

    gen iso3_d2 = iso3_d
    replace iso3_d2 = "EUN" if eu_d == 1
    merge m:1 iso3_o iso3_d2 using "$ROOT/Data/Int/GSP_2012_data"
    drop if _merge == 2
    replace delta = 1 if delta == .
    replace delta_ind = 1 if delta_ind == .
    replace delta_wt = 1 if delta_wt == .
    replace gsp_rate = 1 if gsp_rate == .
    replace mfn_rate = 1 if mfn_rate == .

    keep    iso3_o iso3_d  gsp_d_d  acp_to_eu delta delta_ind delta_wt  gsp_rate mfn_rate
    rename (delta delta_ind delta_wt) (gsp_delta gsp_delta_ind gsp_delta_wt )
    saveold "$ROOT/Data/Int/GSP", replace

    // Origin and Destination Wealth for G
    use "$ROOT/Data/Int/stack_data_B", clear
    keep if year == 2000
    keep rich_dest  rich_orig  iso3_o iso3_d gdpcap_o  
    merge 1:1 iso3_o iso3_d using  "$ROOT/Data/Int/GSP"
    replace gsp_d_d = 0 if gsp_d_d == .
    duplicates drop
    saveold "$ROOT/Data/Int/wealth2000B", replace

    // EORA Data
    insheet using  "$ROOT/Data/EORA/Output/trade_flows_all_sectors.csv", clear

    // Something is wrong with ETHIOPIA!
    replace destination = "ETH" if origin == "ETH" & destination == "ROW"
    replace x = x*10 if origin == "ETH" & destination == "ETH"
    rename x X_ij_EORA
    rename (origin destination) (iso3_o iso3_d)
    gcollapse (sum) X_ij_EORA, by( iso3_d iso3_o year)
    saveold  "$ROOT/Data/Int/trade_flows_all_sectors", replace

    // Goal impute n_ij
    use  "$ROOT/Data/Int/stack_data_B", clear
    keep if year == ${year}

    cap drop _merge
    merge 1:1 iso3_o iso3_d year using "$ROOT/Data/Int/trade_flows_all_sectors"
    drop if _merge == 1
    count

    replace X_ij = X_ij_EORA
    assert X_ij > 0
    replace lX_ij = log(X_ij)

    drop if ldistw == .
    gen lGDP_diff = log(gdpcap_o/gdpcap_d)
    regress lX_ij ldistw fta_wto comcur comrelig  sibling tdiff contig lGDP_diff i.o i.d
    predict  lX_ij_predicted, xb
    gen X_ij_imputed = cond(X_ij!=.,X_ij,exp(lX_ij_predicted))
    gen lX_ij_imputed = log(X_ij_imputed)

    preserve
      collapse (sum) X_ij, by(iso3_o)
      list if inlist(iso3_o,"ROW","TOT")
    restore

    // Drop MMR
    drop if inlist(iso3_o,"ROW","TOT","MMR","BWA")
    drop if inlist(iso3_d,"ROW","TOT","MMR","BWA")

    preserve
      gen have = 1 if ln_ij != .
      collapse (sum) X_ij (count)  count = X_ij, by(have)
      egen X_ij_total = total(X_ij)
      egen count_total = total(count)
      gen X_ij_share = X_ij/X_ij_total
      gen count_share = count/count_total
      list
      outsheet using "$ROOT/Output/Statistics${year}Big.csv", replace comma
      // This is for the share
    restore

    merge 1:1 iso3_o iso3_d using "$ROOT/Data/Int/wealth2000B", nogen keep(master matched)


    // First impute using distance and trade flows
    cap drop ln_ij_predicted
    /* regress ln_ij ldistw lX_ij_imputed i.o i.d */
    regress ln_ij ldistw lX_ij_imputed fta_wto comcur comrelig  sibling tdiff contig lGDP_diff i.o i.d rich_dest#rich_orig
    predict  ln_ij_predicted, xb
    summ ln_ij
    /* replace ln_ij_predicted = `r(max)' if ln_ij_predicted > `r(max)' */
    replace ln_ij_predicted = . if iso3_o == iso3_d
    summ ln_ij_predicted
    gen double n_ij_imputed = cond(n_ij!=. & n_ij!=0,n_ij,exp(ln_ij_predicted))
    gen double ln_ij_imputed = cond(ln_ij!=. & n_ij!=0 ,ln_ij,(ln_ij_predicted))


    summ n_ij if iso3_o == iso3_d, d
    replace n_ij_imputed = `r(p50)' if iso3_o == iso3_d & n_ij == .

    summ ln_ij if iso3_o == iso3_d, d
    replace ln_ij_imputed = `r(p50)' if iso3_o == iso3_d & n_ij == .


    keep iso3_o iso3_d n_ij ln_ij n_ij_imputed ln_ij_imputed X_ij lX_ij X_ij_imputed lX_ij_imputed
    sort iso3_o iso3_d
    summ n_ij_imputed, d
    replace n_ij_imputed = `r(p1)' if n_ij_imputed <= `r(p1)'
    saveold "$ROOT/Data/Int/n_ij_imputed", replace


    use "$ROOT/Data/Int/n_ij_imputed", clear
    drop if inlist(iso3_o,"ROW","TOT","MMR")
    drop if inlist(iso3_d,"ROW","TOT","MMR")
    keep n_ij_imputed iso3_o iso3_d
    reshape wide n_ij_imputed, i(iso3_o) j(iso3_d) string
    keep iso3_o
    outsheet using "$ROOT/Data/Int/WIOD_sampleB/l_i_${year}B.csv", replace non comma

    

    // Matrix: Column is destination; Row is Origin

    // Get x_ij and y_ij
    use "$ROOT/Data/Int/n_ij_imputed", clear
    merge 1:1 iso3_o iso3_d using "$ROOT/Data/Int/wealth2000B", nogen keep(master matched)

    merge m:1 iso3_o using "$ROOT/Data/Int/estimation_sample_${year}", keep(master matched) gen(sample)


    drop X_ij n_ij l*

    bys iso3_d: egen  E_j = total(X_ij)
    gen x_ij = X_ij/E_j

    bys iso3_o: egen  Y_i = total(X_ij)
    gen y_ij = X_ij/Y_i
    gen G_ij = .
    replace G_ij = 0 if rich_dest == 1 & rich_orig == 1
    replace G_ij = 1 if rich_dest == 0 & rich_orig == 1
    replace G_ij = 2 if rich_dest == 1 & rich_orig == 0
    replace G_ij = 3 if rich_dest == 0 & rich_orig == 0
    tab G_ij, missing


    gen G2_ij = 0
    replace G2_ij = 1 if rich_orig == 1
    replace G2_ij = 1 if rich_orig == 1
    tab G2_ij, missing

    gen G4_ij = .
    replace G4_ij = 1 if rich_orig == 1
    replace G4_ij = 0 if rich_orig == 0 
    replace G4_ij = 3 if rich_orig == 0 & gdpcap_o <= 2995 // 7.3140e+03
    replace G4_ij = 2 if rich_orig == 0 & gdpcap_o <= 755 // 2.9217e+03
    tab G4_ij, missing


    preserve
        keep iso3_o sample
        sort iso3_o 
        duplicates drop
        drop iso*
        outsheet using "$ROOT/Data/Int/WIOD_sampleB/Sample_${year}B.csv", replace non nol comma
        count
    restore


    preserve
        keep iso3* gsp_delta_ind
        sort iso3_o iso3_d
        reshape wide gsp_delta_ind, i(iso3_o) j(iso3_d) string
        drop iso*
        outsheet using "$ROOT/Data/Int/WIOD_sampleB/GSP_delta_ind_${year}B.csv", replace non nol comma
        count
    restore

    preserve
        keep iso3* gsp_delta
        sort iso3_o iso3_d
        reshape wide gsp_delta, i(iso3_o) j(iso3_d) string
        drop iso*
        outsheet using "$ROOT/Data/Int/WIOD_sampleB/GSP_delta_${year}B.csv", replace non nol comma
        count
    restore

    preserve
        keep iso3* mfn_rate
        sort iso3_o iso3_d
        reshape wide mfn_rate, i(iso3_o) j(iso3_d) string
        drop iso*
        outsheet using "$ROOT/Data/Int/WIOD_sampleB/mfn_rate_${year}B.csv", replace non nol comma
        count
    restore


    preserve
        keep iso3* gsp_d_d
        sort iso3_o iso3_d
        reshape wide gsp_d_d, i(iso3_o) j(iso3_d) string
        drop iso*
        outsheet using "$ROOT/Data/Int/WIOD_sampleB/GSP_${year}B.csv", replace non nol comma
        count
    restore

    preserve
        keep iso3* n_ij_imputed
        sort iso3_o iso3_d
        reshape wide n_ij, i(iso3_o) j(iso3_d) string
        drop iso*
        outsheet using "$ROOT/Data/Int/WIOD_sampleB/n_ij_${year}B.csv", replace non nol comma
        count
    restore

    preserve
        keep iso3* G_ij
        sort iso3_o iso3_d
        reshape wide G_ij, i(iso3_o) j(iso3_d) string
        drop iso*
        outsheet using "$ROOT/Data/Int/WIOD_sampleB/G_ij_${year}B.csv", replace non nol comma
        count
    restore

    preserve
        keep iso3* G2_ij
        sort iso3_o iso3_d
        reshape wide G2_ij, i(iso3_o) j(iso3_d) string
        drop iso*
        outsheet using "$ROOT/Data/Int/WIOD_sampleB/G2_ij_${year}B.csv", replace non nol comma
        count
    restore


    preserve
        keep iso3* G4_ij
        sort iso3_o iso3_d
        reshape wide G4_ij, i(iso3_o) j(iso3_d) string
        drop iso*
        outsheet using"$ROOT/Data/Int/WIOD_sampleB/G4_ij_${year}B.csv", replace non nol comma
        count
    restore


    preserve
        keep iso3* x_ij
        sort iso3_o iso3_d
        reshape wide x_ij, i(iso3_o) j(iso3_d) string
        drop iso*
        outsheet using "$ROOT/Data/Int/WIOD_sampleB/x_ij_${year}B.csv", replace non nol comma
        count
    restore

    preserve
        keep iso3* X_ij
        sort iso3_o iso3_d
        reshape wide X_ij, i(iso3_o) j(iso3_d) string
        drop iso*
        outsheet using "$ROOT/Data/Int/WIOD_sampleB/XX_ij_${year}B.csv", replace non nol comma
    restore

    preserve
        keep iso3* X_ij_imputed
        sort iso3_o iso3_d
        reshape wide X_ij_imputed, i(iso3_o) j(iso3_d) string
        drop iso*
        outsheet using "$ROOT/Data/Int/WIOD_sampleB/XX_imputed_ij_${year}B.csv", replace non nol comma
    restore

    preserve
        keep iso3* y_ij
        sort iso3_o iso3_d
        reshape wide y_ij, i(iso3_o) j(iso3_d) string
        drop iso*
        outsheet using "$ROOT/Data/Int/WIOD_sampleB/y_ij_${year}B.csv", replace non nol comma
    restore

    preserve
        keep iso3_d E_j
        rename iso3_d iso3
        rename E_j E_i
        duplicates drop
        tempfile T
        saveold `T', replace
    restore

    keep iso3_o Y_i
    rename iso3 iso3
    duplicates drop
    merge 1:1 iso3  using `T', nogen
    sort iso
    gen kappa_i = Y_i/E_i
    drop iso*
    outsheet using "$ROOT/Data/Int/WIOD_sampleB/balance_i_${year}B.csv", replace non nol comma


    use "$ROOT/Data/Int/stack_data_B", clear
    keep if year == 2012
    /* keep if iso3_d == iso3_o */
    tempfile A 
    saveold `A', replace

    use "$ROOT/Data/Int/n_ij_imputed", clear
    /* keep if iso3_d == iso3_o */
    merge 1:1 iso3_d iso3_o using `A',  keepusing(gdp_o gdpcap_o area_o lN_ii)   keep(master match) 
    assert _merge == 3

    /* use $ROOT/Data/Int/stack_data_${year}, clear */

    // gen lN_ii = log(N_ii)
    gen lgdp_o = log(gdp_o)
    gen lgdpcap_o = log(gdpcap_o)
    gen larea_o   = log(area_o)

    regress lN_ii   lgdp_o lgdpcap_o larea_o gdp_o gdpcap_o area_o
    predict lN_ii_imputed, xb
    gen N_ii_imputed = exp(lN_ii_imputed)
    order N_ii N_ii_imputed
    // scatter N_ii_imputed N_ii

    keep iso3* N_ii_imputed
    format N* %15.0fc
    replace N_ii_imputed = N_ii if N_ii != .
    tempfile N_ii_imputed
    saveold `N_ii_imputed', replace
    sort iso3_o iso3_d
    summ N_ii_imputed, d
    replace N_ii_imputed = `r(p1)' if N_ii_imputed <= `r(p1)'
    reshape wide N_ii_imputed, i(iso3_o) j(iso3_d) string
    outsheet using "$ROOT/Data/Int/WIOD_sampleB/N_ii_${year}B_lab.csv", replace non comma
    preserve
        drop iso*
        /* gcollapse (sum) * */
        outsheet using "$ROOT/Data/Int/WIOD_sampleB/N_ii_${year}B.csv", replace non nol comma
    restore
    keep iso3_o
    gen i = _n
    outsheet using "$ROOT/Data/Int/WIOD_sampleB/Names_${year}B.csv", replace non nol comma




end

// Run Routine:
MAIN

// EOF
