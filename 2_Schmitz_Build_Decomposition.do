/* Created by RM on 2018.05.04 to calculate TFP using Hsieh Klenow measures */
*TRIMMING MORE


set more off
clear

global raw "/Users/russellmorton/Desktop/SRA 2018/STATA/Data/Raw"

global temp "/Users/russellmorton/Desktop/SRA 2018/STATA/Data/Temp"

global cleaned "/Users/russellmorton/Desktop/SRA 2018/STATA/Data/Cleaned"

global out "/Users/russellmorton/Desktop/SRA 2018/STATA/Output/Schmitz Productivity Decompositions"

use "$raw/firms", clear

*Check if there are any missing years

drop if routput == . | routput == 0 | wages == . | wages == 0 | capn == . | capn == 0 
drop if rmata == . | rmata == 0

xtset firm year

bys firm: egen min_firm_year = min(year)

bys firm: egen max_firm_year = max(year)

bys firm: egen count_firm_years = count(year)

drop if count_firm_years == 1

g max_less_min_year = (max_firm_year - min_firm_year) + 1

g diff_years = max_less_min_year - count_firm_years

drop if diff_year > 0

g routput_1996 = 0
replace routput_1996 = routput if year == 1996
g wages_1996 = 0
replace wages_1996 = wages if year == 1996
g capn_1996 = 0
replace capn_1996 = capn if year == 1996
g rmata_1996 = 0
replace rmata_1996 = rmata if year == 1996

bys sector: egen sum_routput_1996 = sum(routput_1996)
bys sector: egen sum_wages_1996 = sum(wages_1996)
bys sector: egen sum_capn_1996 = sum(capn_1996)
bys sector: egen sum_rmata_1996 = sum(rmata_1996)

g pre_theta_labor_1996 = sum_wages_1996 / sum_routput_1996
g pre_theta_cap_1996  =  sum_capn_1996 / sum_routput_1996
g pre_theta_materials_1996 = sum_rmata_1996 / sum_routput_1996

g theta_labor_1996 =  pre_theta_labor_1996 / (pre_theta_labor_1996 + pre_theta_cap_1996 + pre_theta_materials_1996)
g theta_cap_1996 =  pre_theta_cap_1996 / (pre_theta_labor_1996 + pre_theta_cap_1996 + pre_theta_materials_1996)
g theta_mat_1996 =  pre_theta_materials_1996 / (pre_theta_labor_1996 + pre_theta_cap_1996 + pre_theta_materials_1996)

*Now Compute Sectoral Level Variables

bys year sector: egen sectoral_routput = sum(routput)
bys year sector: egen sectoral_wages = sum(wages)
bys year sector: egen sectoral_capn = sum(capn)
bys year sector: egen sectoral_mat = sum(rmata)

g sectoral_tfp = sectoral_routput / (sectoral_wages ^ theta_labor_1996 * sectoral_capn ^ theta_cap_1996 * sectoral_mat ^ theta_mat_1996 )
g sectoral_labor_prod = sectoral_routput / (sectoral_wages ^ theta_labor_1996 )
g sectoral_cap_prod = sectoral_routput / (sectoral_capn ^ theta_cap_1996 )
g sectoral_mat_prod = sectoral_routput / (sectoral_mat ^ theta_mat_1996)

egen tag_for_graph = tag(year sector)

*Now standardize to the first year in the sample

save "$temp/Balanced Panel Only", replace

*preserve

keep if tag_for_graph == 1
xtset sector year

bys sector: egen min_year = min(year)

foreach var in tfp labor_prod cap_prod mat_prod {
	g pegged_`var' = 1 if year == min_year
	forv i = 1(1)20 {
		replace pegged_`var' =  L1.pegged_`var' * (1 + ( (sectoral_`var' - L1.sectoral_`var') /  L1.sectoral_`var') ) if pegged_`var' == . & L1.pegged_`var' != . & sector == L1.sector
	}
}	
	
label var pegged_tfp "Pegged Sectoral TFP"
label var pegged_labor_prod "Pegged Sectoral Labor Productivity"
label var pegged_cap_prod "Pegged Sectoral Capital Productivity"
label var pegged_mat_prod "Pegged Sectoral Materials Productivity"

egen min_sector = min(sector)
egen max_sector = max(sector)
local min_sector = min_sector
local max_sector = max_sector

forv s = `min_sector'(1)`max_sector' {
*local s = 7

	di "enters sector loop"
	preserve
	keep if sector == `s'
	local secname = secname	
	
	capture drop min_year max_year
	egen sector_min_year = min(year)
	egen sector_max_year = max(year)
	local minyear = sector_min_year
	local maxyear = sector_max_year	
	
	di "makes it to graph"
	
	di "`minyear' is minyear and `maxyear' is maxyear"
		
twoway (line pegged_tfp year) ///
	   (line pegged_labor_prod year) ///
	   (line pegged_cap_prod year) ///
	   (line pegged_mat_prod year, ///
	   xmtick(`minyear'(1)`maxyear') xlab(`minyear'(1)`maxyear') ///
        xtitle("Year") ytitle("Productivity (Pegged to `minyear' = 1)", size(medsmall)) title("Productivity over Time for Ghana `secname' Sector")  ///
		note("Source: CSAE Ghana RPED/GMES Data") legend(size(vsmall)) )

	graph export "$out/Productivty Over Time for `secname'.pdf", as (pdf) replace
		
	restore	
		
} 
	   

*Do Decomposition

use "$temp/Balanced Panel Only", clear

bys sector: egen sector_min_year = min(year)
bys sector: egen sector_max_year = max(year)
	
g first_year = 0
g continuing = 0
g entrant = 0
g exit = 0 
g last_year = 0

forv s = `min_sector'(1)`max_sector' {
	egen pre_minyear = max(sector_min_year) if sector == `s'
	egen minyear = max(pre_minyear)
	egen pre_maxyear = max(sector_max_year) if sector == `s'
	egen maxyear = max(pre_maxyear)
 
	local minyear = minyear
	local maxyear = maxyear	
	
	drop pre_minyear pre_maxyear minyear maxyear

	forv y = `minyear'(1)`maxyear' {
		replace first_year = 1 if year == `minyear' & sector == `s'
		replace continuing = 1 if year > min_firm_year & sector == `s' & year < max_firm_year
		replace entrant = 1 if year != `minyear' & sector == `s' & year == min_firm_year
		replace exit = 1 if year != `maxyear' & sector == `s' & year == max_firm_year & year != `minyear'
		replace last_year = 1 if year == `maxyear'
		}
	}


g check_sum = first_year + continuing + entrant + exit + last_year
bys firm: egen max_check_sum = max(check_sum)
	
sort sector firm year
*browse sector firm year first_year continuing entrant exit last_year check_sum if max_check_sum > 1

*Now Make Table II

save "$temp/pre_table_II", replace

*Now create shell for merging
use "$temp/pre_table_II", clear

xtset firm year
g lag_year = L1.year
keep if year == lag_year + 1 | year == sector_min_year
bys sector year: egen min_firm_for_keep = min(firm)
keep if min_firm_for_keep == firm

keep sector year sector_min_year sector_max_year
sort sector year

keep sector year
save "$temp/sector year shell for merge", replace

use "$temp/pre_table_II", clear

keep continuing exit entrant routput wages firm sector year sector_min_year sector_max_year

save "$temp/data for merge", replace

use "$temp/pre_table_II", clear

keep sector firm
duplicates drop

save "$temp/sector firm shell for merge", replace

use "$temp/sector year shell for merge", clear

joinby sector using "$temp/sector firm shell for merge"

save "$temp/shell for merge", replace

use "$temp/shell for merge", clear

merge 1:1 sector firm year using "$temp/data for merge"

bys year sector: egen sum_year_sector_routput = sum(routput)
bys year sector: egen sum_year_sector_wages = sum(wages)
bys year sector: egen count_continuing = sum(continuing)
bys year sector: egen count_firms_sector_year = count(firm)

g s_it = wages / sum_year_sector_wages
g pi_it = routput / wages 
replace s_it = 0 if s_it == .
replace pi_it = 0 if pi_it == .

bys year sector: g pi = sum_year_sector_routput / sum_year_sector_wages
xtset firm year

g overall_growth = D1.pi
g delta_pi_it = D1.pi_it
g delta_s_it = D1.s_it

xtset firm year
g pre_within_mine = 0
replace pre_within_mine = L1.s_it * delta_pi_it if continuing == 1
bys year sector: egen within_mine = sum(pre_within_mine)

xtset firm year
g pre_between_mine = 0
replace pre_between_mine = (L1.pi_it - L1.pi) * delta_s_it if continuing == 1
bys year sector: egen between_mine = sum(pre_between_mine)

xtset firm year
g pre_cross_mine = 0
replace pre_cross_mine = delta_pi_it * delta_s_it if continuing == 1
bys year sector: egen cross_mine = sum(pre_cross_mine)

xtset firm year
g pre_exit_prod = 0
replace pre_exit_prod = L1.s_it * (L1.pi_it - L1.pi) if exit == 1
bys year sector: egen exit_prod = sum(pre_exit_prod)

xtset firm year
g pre_entrant_prod = 0
replace pre_entrant_prod = s_it * (pi_it - L1.pi) if entrant == 1
bys year sector: egen entrant_prod = sum(pre_entrant_prod)

*keep if sector == 2 & year <= 1992
*export excel using "$temp/scratch.xlsx", replace first(var)
*1994sector 1

keep if sector == 1 & year <= 1994 & year >= 1993
export excel using "$temp/debug sector 1 1994.xlsx", replace first(var)


*preserve
xtset firm year
g lag_year = L1.year
keep if year == lag_year + 1 | year == sector_min_year
bys sector year: egen min_firm_for_keep = min(firm)
keep if min_firm_for_keep == firm

sort sector year
keep sector year overall_growth within_mine between_mine cross_mine exit_prod entrant_prod

g check_sum_prod =  within_mine + between_mine + cross_mine ///
						- exit_prod + entrant_prod


sort sector year
browse sector year overall_growth check_sum_prod		

*keep if sector == 1 & year <= 1993
