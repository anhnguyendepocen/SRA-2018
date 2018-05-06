/* Created by RM on 2018.05.04 to calculate TFP using Hsieh Klenow measures */
*TRIMMING MORE


set more off
clear

global raw "/Users/russellmorton/Desktop/SRA 2018/STATA/Data/Raw"

global temp "/Users/russellmorton/Desktop/SRA 2018/STATA/Data/Temp"

global cleaned "/Users/russellmorton/Desktop/SRA 2018/STATA/Data/Cleaned"

global out "/Users/russellmorton/Desktop/SRA 2018/STATA/Output/Hsieh Klenow TFP Measures"

use "$raw/firms", clear

local sigma = 3


drop if capn == . | wages == . | routput == .
drop if routput < 1

bys year: egen count_year = count(year)


**Calculate labor share for alpha by sector
*Based on counts above, use 1996

g wages_1996 = 0
replace wages_1996 = wages if year == 1996
g Y_1996 = 0
replace Y_1996 = routput if year == 1996

bys sector: egen sum_wages_1996 = sum(wages_1996)
bys sector: egen sum_Y_1996 = sum(Y_1996)
g alpha_1996 = sum_wages_1996 / sum_Y_1996

**Calculate Marginal Product of Labor and of Capital for Use in TFPR bar

foreach var in capn wages routput {

bys sector year: egen sum_`var'_raw = sum(`var')

}

g tfpr_bar_raw = sum_routput_raw / ( (sum_wages_raw^alpha_1996) * (sum_capn_raw^(1-alpha_1996)) )

g tfpr_i = routput / ( (wages^alpha_1996) * (capn^(1-alpha_1996)) )

g log_tfpr_raw = log(tfpr_i/tfpr_bar_raw)

g tfpq_i= (routput^((`sigma')/(`sigma'-1)))/((capn^alpha_1996) * (wages ^ (1-alpha_1996)))

bys year sector: egen tfpq_bar_raw = mean(tfpq_i)

g log_tfpq_raw = log(tfpq_i/ tfpq_bar_raw)

g drop_flag = 0

foreach prod in tfpr tfpq {
			di "enters prod loop.  prod is `prod'"
			foreach i in 8 92 {
				di "enters values loop.  prod is `prod' and i is `i'"
			
				bys year sector: egen pct_log_`prod'_`i'_raw = pctile(log_`prod'_raw), p(`i')
				
				di "generates pctile"
				
			if `i' < 10 {
				replace drop_flag = 1 if log_`prod'_raw <= pct_log_`prod'_`i'_raw
				}
			else {
				replace drop_flag = 1 if log_`prod'_raw >= pct_log_`prod'_`i'_raw
				}
			}
		}

		
save "$temp/pre drop", replace

use "$temp/pre drop", clear

drop if drop_flag == 1

foreach var in capn wages routput {

bys sector year: egen sum_`var' = sum(`var')

}

g tfpr_bar = sum_routput / ( (sum_wages^alpha_1996) * (sum_capn^(1-alpha_1996)) )
g log_tfpr = log(tfpr_i/tfpr_bar)


bys year sector: egen tfpq_bar = mean(tfpq_i)
g log_tfpq = log(tfpq_i/ tfpq_bar)

su log_tfpr, d
su log_tfpq, d

foreach prod in tfpr tfpq {
	bys year: egen sd_log_`prod' = sd(log_`prod')
		foreach j in 10 25 75 90 {
			bys year: egen pct_`j'_log_`prod' = pctile(log_`prod'), p(`j')
		}
	g disp_75_25_log_`prod' = pct_75_log_`prod' - pct_25_log_`prod'
	g disp_90_10_log_`prod' = pct_90_log_`prod' - pct_10_log_`prod'
	
	bys year: egen nobs_`prod' = count(log_`prod')
	}

	 *kdensity log_tfpr if year == 1996
	 

save "$temp/All_Years_Log_TFPR_TFPQ", replace

	 
* pick 1992 1996 and 1999
keep if inlist(year,1992,1996,1999)


///*******
*EFFICIENT OUTPUT : ECONOMY WIDE
********///


*Identify sectoral shares

g routput_1996 = 0
replace routput_1996 = routput if year == 1996
bys sector: egen sectoral_output = sum(routput_1996)

egen total_output_1996 = sum(routput_1996)

g theta_1996 = sectoral_output / total_output_1996



*Calculate efficient output
/*OLD

bys year sector: egen tfpq_effic = sum(tfpq_i^(`sigma'-1))
bys year sector: egen tfpq_effic_count = count(tfpq_i)
g tfpq_bar_effic = ((tfpq_effic)/tfpq_effic_count)^(1/(`sigma'-1))
g sectoral_effic_Y = tfpq_bar_effic * (capn^(1-alpha_1996)) * (wages ^ alpha_1996)
*g sectoral_effic_Y = tfpq_bar * (capn^(1-alpha_1996)) * (wages ^ alpha_1996)
g weighted_sectoral_Y = (sectoral_effic_Y / 100000000000) ^ theta_1996 
g ln_weighted_sectoral_Y = ln(weighted_sectoral_Y)
bys year: egen ln_economy_efficient_Y = sum(ln_weighted_sectoral_Y)
g efficient_Y = exp(ln_economy_efficient_Y)

*Use HK equation 15 for sectoral TFP
g actual_Y_i = tfpq_i * (capn^(1-alpha_1996)) * (wages ^ alpha_1996)
g weighted_actual_Y = (actual_Y / 100000000000) ^ theta_1996
g ln_weighted_actual_Y = ln(weighted_actual_Y)
bys year: egen ln_economy_actual_Y = sum(ln_weighted_actual_Y)
g actual_Y = exp(ln_economy_actual_Y)

g percent_of_effic = actual_Y/efficient_Y
*/
*Efficient
bys year sector: egen tfpq_effic = sum(tfpq_i^(`sigma'-1))
bys year sector: egen tfpq_effic_count = count(tfpq_i)
*g tfpq_bar_effic = ((tfpq_effic)/tfpq_effic_count)^(1/(`sigma'-1))
g tfpq_bar_effic = ((tfpq_effic))^(1/(`sigma'-1))
g effic_sector_output_theta = (tfpq_bar_effic * (sum_capn^(1-alpha_1996)) * (sum_wages ^ alpha_1996) )^theta_1996

*Actual
bys year sector: egen pre_tfp_s_actual = sum((tfpq_i*tfpr_bar / tfpr_i)^(`sigma'-1))
bys year sector: egen pre_tfp_s_actual_count = count(tfpq_i)
g tfp_s = ((pre_tfp_s_actual)/pre_tfp_s_actual_count)^(1/(`sigma'-1))
g actual_sector_output_theta = (tfp_s * (sum_capn^(1-alpha_1996)) * (sum_wages ^ alpha_1996) )^theta_1996

*Now Sum
egen tag_for_output = tag(sector year)

foreach type in effic actual {
g ln_`type'_output_sector = 0
replace ln_`type'_output_sector = ln(`type'_sector_output_theta) if tag == 1
bys year: egen ln_economy_`type'_output = sum(ln_`type'_output_sector)
g economy_`type'_output = exp(ln_economy_`type'_output)
}

tab year ln_economy_effic_output
tab year ln_economy_actual_output

g output_gains_from_tfp = 100*((economy_effic_output/economy_actual_output)-1)
tab year output_gains_from_tfp


///*** MAKE OUTPUT ***///

g tfpq_dispersion_graph = log(tfpq_i * (tfpq_effic_count^(1/(`sigma'-1))) / tfpq_bar_effic)
kdensity tfpq_dispersion_graph if year == 1996



