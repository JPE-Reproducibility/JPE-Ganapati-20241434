// HS2 Level

// Get Rich Countries
capture program drop setupdata_Rich
program define setupdata_Rich
  use "$ROOT/Data/Int/CEPII_Grav_temp" , clear
  gen gdpcap_d_2000_ = gdpcap_d if year == 2000
  gen gdpcap_o_2000_ = gdpcap_o if year == 2000
  bys iso3_d: egen gdpcap_d_2000 = max(gdpcap_d_2000_)
  bys iso3_o: egen gdpcap_o_2000 = max(gdpcap_o_2000_)
  drop gdpcap_d_2000_ gdpcap_o_2000_
  gen rich_dest = cond(gdpcap_d_2000>9000,1,0)
  gen rich_orig = cond(gdpcap_o_2000>9000,1,0)
  keep iso3_o iso3_d  rich_dest rich_orig gdpcap_d_2000 gdpcap_o_2000
  duplicates drop
  saveold "$ROOT/Data/Int/CEPII_rich", replace
end

// Setup HS 1-digit
capture program drop setupdata_HS1
program define setupdata_HS1
  // EDD (Need an impute step here)
  use "$ROOT/Data/EDD_Data/CYH2D_manuf/CYH2D_manuf.dta", clear
  keep c h2 d y A1 A6i A6ii
  keep if inlist(y,2010,2012)
  rename c iso3_o
  rename d iso3_d
  rename y year
  rename A1 N_ij_h2
  rename A6i x_ij_h2

  gen hs2 = substr(h2,1,2)
  gen X_ij_h2 = N_ij_h2*x_ij_h2

  gen w = N_ij_h2

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

  gcollapse (rawsum) N_ij_h2 (mean) x_ij_h2 [aw=N_ij_h2], by(iso3_o iso3_d year hs1)

  egen group_oy = group(iso3_o year h)
  egen group_dy = group(iso3_d year h)
  egen group_od = group(iso3_d year iso3_o)

  // Merge in Tariff Data
  tab h
  gen h = hs1
  merge m:1 iso3_o iso3_d hs1 year using "$ROOT/Data/TRAINS/collapsed/T_all_collapse_hs1.dta", keep(master matched)
  drop _merge

  merge m:1 iso3_o iso3_d hs1 year using "$ROOT/Data/TRAINS/collapsed/Teti2012_all_collapse_hs1.dta", keep(master matched)
  drop _merge

  replace simpleAHS_w = 0 if iso3_o == iso3_d
  replace simpleAHS_uw = 0 if iso3_o == iso3_d
  replace TsimpleAHS_uw = 0 if iso3_o == iso3_d

  // Merge in CEPII Gravity
  merge m:1 iso3_o iso3_d year using "$ROOT/Data/Int/CEPII_Grav_temp", keep(master matched)
  drop _merge
  merge m:1 iso3_o iso3_d using "$ROOT/Data/Int/CEPII_rich", nogen

  // Setup Variabes
  gen lx_ij_h1 = log(x_ij_h2)
  gen ldist = log(dist)
  gen ldistw = log(distw)
  gen lN_ij_h1 = log(N_ij_h2)
  gen lX_ij_h1 = lx_ij_h1 + lN_ij_h1
  global year = 2012
  keep if year == $year
  destring h, replace

  // The next lines is a rough approximation
  // We need a factor to scale the data
  // Simply use 10% here
  bys iso3_o : egen max_N_ij =max(N_ij_)
  gen ln_ij = log(N_ij/(max_N_ij*10))

  // Clean Up Variables for Output
  tab contig, missing
  replace contig = 1 if iso3_o == iso3_d
  tab comlang_ethno, missing
  replace fta_wto = 0 if fta_wto == .
  replace fta_wto = 1 if iso3_o == iso3_d
  gen dist_own  = distw if  iso3_o == iso3_d
  bys iso3_o: egen dist_ii = min(dist_own)
  bys iso3_d: egen dist_jj = min(dist_own)
  drop dist_own

  gen tariff_use = simpleAHS_w/100
  gen ltariff_use = log(1+tariff_use)
  gen ltariff_use2 = ltariff_use

  gen tariff_useUW = simpleAHS_uw/100
  gen ltariff_useUW = log(1+tariff_useUW)
  gen ltariff_use2UW = ltariff_useUW

  gen Ttariff_useUW = TsimpleAHS_uw/100
  gen Tltariff_useUW = log(1+Ttariff_useUW)
  gen Tltariff_use2UW = Tltariff_useUW

  gen ldist_use = ldistw
  gen lx_bar = lx_ij_h

  gen dummy1 = 1
  gen dummy2 = 1
  gen dummy3 = ln(pop_o )*ln(pop_d)
  gen dummy4 = log(dist_ii)*log(dist_jj)
  preserve
    keep iso3_o
    duplicates drop
    gen sectoral = "1"
    saveold  "$INT/sectoral_sample", replace
    count
  restore

end

// Subrountine to setup 1-digit
capture program drop make_HS1
program define make_HS1
syntax, keeph(integer)
global keeph = `keeph'

  // Unweighted Tariffs
  preserve
    if $keeph != 0 {
      keep if h == $keeph
    } 
    else {
        keep if inlist(h,3,4,5,6,7,8,9,10,11,12,13,14)
        /* 5, 13 is super wide */
    }
    eststo A1n: reghdfe lx_ij_h i.rich_orig#(i.rich_dest) #(c.ldistw c.ltariff_use2UW)#i.h , absorb(group_oy group_dy h) cluster(group_oy group_dy ) 
    keep if e(sample)


    tostring group_oy, gen(tgroup_oy)
    tostring group_dy, gen(tgroup_dy)
    encode tgroup_oy, gen(oN)
    encode tgroup_dy, gen(dN)

    assert  ldistw != .
    assert  lN_ij != .
    keep if ltariff_use2UW != .

    keep  oN dN lN_ij lx_bar ltariff_useUW gdpcap_o gdpcap_d  rich_orig rich_dest ln_ij gdp_o gdp_d colony contig fta_wto comcur comlang_ethno  dummy1 dummy2 dummy3 dummy4 ldist_use ltariff_use2UW h
    order oN dN lN_ij lx_bar ltariff_useUW gdpcap_o gdpcap_d  rich_orig rich_dest ln_ij gdp_o gdp_d colony contig fta_wto comcur comlang_ethno  dummy1 dummy2 dummy3 dummy4 ldist_use ltariff_use2UW h
    sort h
    outsheet using "$INT/WIOD_sampleB/${year}_Tuw_h${keeph}.csv", replace non nol comma
    count
    tab rich_orig rich_dest
  restore


  // Unweighted Tariffs - Teti
  preserve
    if $keeph != 0 {
      keep if h == $keeph
    } 
    else {
        keep if inlist(h,3,4,5,6,7,8,9,10,11,12,13,14)
    }
    eststo A1n: reghdfe lx_ij_h i.rich_orig#(i.rich_dest) #(c.ldistw c.Tltariff_use2UW)#i.h , absorb(group_oy group_dy h) cluster(group_oy group_dy ) 
    keep if e(sample)

    tostring group_oy, gen(tgroup_oy)
    tostring group_dy, gen(tgroup_dy)
    encode tgroup_oy, gen(oN)
    encode tgroup_dy, gen(dN)

    assert  ldistw != .
    assert  lN_ij != .
    keep if Tltariff_use2UW != .

    keep  oN dN lN_ij lx_bar Tltariff_useUW gdpcap_o gdpcap_d  rich_orig rich_dest ln_ij gdp_o gdp_d colony contig fta_wto comcur comlang_ethno  dummy1 dummy2 dummy3 dummy4 ldist_use Tltariff_use2UW h
    order oN dN lN_ij lx_bar Tltariff_useUW gdpcap_o gdpcap_d  rich_orig rich_dest ln_ij gdp_o gdp_d colony contig fta_wto comcur comlang_ethno  dummy1 dummy2 dummy3 dummy4 ldist_use Tltariff_use2UW h
    sort h
    outsheet using "$INT/WIOD_sampleB/${year}_Tetiuw_h${keeph}.csv", replace non nol comma
    count
    tab rich_orig rich_dest
  restore

end



// Routine
setupdata_Rich
setupdata_HS1
make_HS1, keeph(0) 
 forvalues hh = 3(1)17 {
 cap make_HS1, keeph(`hh')
}


