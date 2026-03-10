local path "$ROOT/Data/China_Exports"

*********************************************************************************************************************************************************************
*********************************************************************************************************************************************************************
*********************************************************************************************************************************************************************
*********************************************************************************************************************************************************************
****YOU SHOULD BE ABLE TO RUN STARTING HERE
***Checking that all countries are covered.
use "`path'/Raw/exports_2000.dta",clear
append using "`path'/Raw/exports_2001.dta"
append using "`path'/Raw/exports_2002.dta"
append using "`path'/Raw/exports_2003.dta"
append using "`path'/Raw/exports_2004.dta"
append using "`path'/Raw/exports_2005.dta"
append using "`path'/Raw/exports_2006.dta"
merge m:1 fpart using "`path'/Raw/fpart_concordance.dta"
tab fpart if _merge==1
keep if _merge==3
keep year n value isocode d
save "`path'/Final/exports_2000_2006.dta",replace


use "`path'/Raw/exports_2007.dta",clear
append using "`path'/Raw/exports_2008.dta"
append using "`path'/Raw/exports_2009.dta"
merge m:1 orig using "`path'/Raw/fpart_concordance.dta"
tab orig if _merge==1
keep if _merge==3
keep year n value isocode d
save "`path'/Final/exports_2007_2010.dta",replace

use "`path'/Raw/exports_2010.dta",clear
append using "`path'/Raw/exports_2011.dta"
append using "`path'/Raw/exports_2012.dta"
merge m:1 d using "`path'/Raw/fpart_concordance.dta"
tab d if _merge==1
keep if _merge==3



keep year n value isocode d
append using "`path'/Final/exports_2000_2006.dta"
append using "`path'/Final/exports_2007_2010.dta"
drop if isocode==""
sort isocode year
rename d country_name
order year isocode country_name value n
save "`path'/Final/Final Data.dta",replace

