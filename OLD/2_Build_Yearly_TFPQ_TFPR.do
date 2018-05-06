
/* Created by RM on 2018.05.03 to calculate TFP using Hsieh Klenow measures */

set more off
clear

global raw "/Users/russellmorton/Desktop/SRA 2018/STATA/Data/Raw"

global temp "/Users/russellmorton/Desktop/SRA 2018/STATA/Data/Temp"

global cleaned "/Users/russellmorton/Desktop/SRA 2018/STATA/Data/Cleaned"

global out "/Users/russellmorton/Desktop/SRA 2018/STATA/Output/Hsieh Klenow TFP Measures"

use "$raw/firms", clear

/*******
* Calculate TFPR and TFPQ
Question: Use all obs, count firms multiple time?
Use last obs per firm?
Use year with most obs? [MOST CONSERVATIVE]
*******/

*First estimate the capital and labor shares by industry for the Cobb Douglas production fct

drop if capn == . | worker == . | routput == .
drop if routput < 1

bys year: egen count_year = count(year)

*Use 1996 to estimate shares

g cobb_douglas_temp = 0
g cobb_douglas_best = 0
*g cobb_douglas_temp_sse = - 10
g cobb_douglas_best_sse = -10

*levelsof(sector), local(sectors)

forv alpha = 0(.0005)1 {
	replace cobb_douglas_temp = `alpha'
	g outputhat = (capn^`alpha') * (worker ^ (1-`alpha'))
	g squared_resid = 0
	replace squared_resid = (outputhat - routput) ^ 2 if year == 1996 
	bys sector: egen sse = sum(squared_resid)
	
	*replace cobb_douglas_temp_sse = sse
	replace cobb_douglas_best_sse = sse if cobb_douglas_best_sse < 0
	replace cobb_douglas_best_sse = sse if sse != . & cobb_douglas_best_sse > sse
	replace cobb_douglas_best = cobb_douglas_temp if cobb_douglas_best_sse == sse
	
	capture drop outputhat squared_resid sse
	}
	
*Use Alt as less biased by survivorship

g alpha_use = cobb_douglas_best

g tfpr_raw = routput /  ((capn^alpha_use) * (worker ^ (1-alpha_use)))
g tfpq_raw = (routput^(3/2))/((capn^alpha_use) * (worker ^ (1-alpha_use)))

*Trim (see pg. 1416)

save "$temp/All Years", replace
use "$temp/All Years", clear

bys year sector: egen sum_tfpr_raw =  sum(tfpr_raw) 
bys year sector: egen count_tfpr_raw =  count(tfpr_raw) 
g avg_sectoral_tfpr_raw = sum_tfpr_raw / count_tfpr_raw

g log_deviation_sectoral_tfpr_raw = log(tfpr_raw/avg_sectoral_tfpr_raw)
bys year sector: egen log_dev_sectoral_tfpr_raw_1pct = pctile(log_deviation_sectoral_tfpr_raw), p(1) 
bys year sector: egen log_dev_sectoral_tfpr_raw_99pct = pctile(log_deviation_sectoral_tfpr_raw), p(99)

bys year sector: egen sum_tfpq_raw =  sum(tfpq_raw) 
bys year sector: egen count_tfpq_raw =  count(tfpq_raw) 
g avg_sectoral_tfpq_raw = sum_tfpq_raw / count_tfpq_raw

g log_deviation_sectoral_tfpq_raw = log(tfpq_raw/avg_sectoral_tfpq_raw)
bys year sector: egen log_dev_sectoral_tfpq_raw_1pct = pctile(log_deviation_sectoral_tfpq_raw), p(1)
bys year sector: egen log_dev_sectoral_tfpq_raw_99pct = pctile(log_deviation_sectoral_tfpq_raw), p(99)

bys year sector: egen count_firms = count(log_deviation_sectoral_tfpq_raw)

g drop_flag =  0
replace drop_flag = 1 if log_deviation_sectoral_tfpr_raw <= log_dev_sectoral_tfpr_raw_1pct 
replace drop_flag = 1 if log_deviation_sectoral_tfpr_raw >= log_dev_sectoral_tfpr_raw_99pct 
replace drop_flag = 1 if log_deviation_sectoral_tfpq_raw <= log_dev_sectoral_tfpq_raw_1pct
replace drop_flag = 1 if log_deviation_sectoral_tfpq_raw >= log_dev_sectoral_tfpq_raw_99pct 

save "$temp/Pre Drop", replace
drop if drop_flag == 1

bys year sector: egen sum_tfpr =  sum(tfpr_raw) 
bys year sector: egen count_tfpr =  count(tfpr_raw) 
g avg_sectoral_tfpr = sum_tfpr/ count_tfpr

bys year sector: egen sum_tfpq =  sum(tfpq_raw) 
bys year sector: egen count_tfpq =  count(tfpq_raw) 
g avg_sectoral_tfpq = sum_tfpq/ count_tfpq

save "$cleaned/firms_TFPR_TFPQ_year", replace

