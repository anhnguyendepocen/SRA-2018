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

**Calculate Marginal Product of Labor and of Capital for Dropping Observations

foreach var in capn wages routput {

bys sector year: egen sum_`var'_raw = sum(`var')

}

g tfpr_bar_raw = sum_routput_raw / ( (sum_wages_raw^alpha_1996) * (sum_capn_raw^(1-alpha_1996)) )

g tfpr_i_raw = routput / ( (wages^alpha_1996) * (capn^(1-alpha_1996)) )

g log_tfpr_raw = log(tfpr_i_raw/tfpr_bar_raw)

g tfpq_bar_raw = (sum_routput_raw^((`sigma')/(`sigma'-1))) / ( (sum_wages_raw^alpha_1996) * (sum_capn_raw^(1-alpha_1996)) )

g tfpq_i_raw = (routput^((`sigma')/(`sigma'-1)))/((capn^alpha_1996) * (wages ^ (1-alpha_1996)))

g log_tfpq_raw = log(tfpq_i_raw/tfpq_bar_raw)

*bys year sector: egen tfpq_bar_raw = mean(tfpq_i)

g drop_flag_tfpr = 0
g drop_flag_tfpq = 0


foreach prod in tfpr tfpq {
			di "enters prod loop.  prod is `prod'"
			foreach i in 8 92 {
				di "enters values loop.  prod is `prod' and i is `i'"
			
				bys year sector: egen pct_log_`prod'_`i'_raw = pctile(log_`prod'_raw), p(`i')
				
				di "generates pctile"
				
			if `i' < 10 {
				replace drop_flag_`prod' = 1 if log_`prod'_raw <= pct_log_`prod'_`i'_raw
				}
			else {
				replace drop_flag_`prod'= 1 if log_`prod'_raw >= pct_log_`prod'_`i'_raw
				}
			}
		}

g drop_flag_comb = drop_flag_tfpr > 0 | drop_flag_tfpq > 0		
		
save "$temp/pre drop", replace

use "$temp/pre drop", clear

drop if drop_flag_comb == 1


***Now Calculate TFPR Measures

g count = 1

foreach var in capn wages routput count {

bys sector year: egen sum_`var' = sum(`var')

}

*First, recalculate alphas

drop wages_1996 Y_1996 sum_wages_1996 sum_Y_1996 alpha_1996 

g wages_1996 = 0
replace wages_1996 = wages if year == 1996
g Y_1996 = 0
replace Y_1996 = routput if year == 1996

bys sector: egen sum_wages_1996 = sum(wages_1996)
bys sector: egen sum_Y_1996 = sum(Y_1996)
g alpha_1996 = sum_wages_1996 / sum_Y_1996


*Now recalculate TFPQ Measures

g tfpq_i = (routput^((`sigma')/(`sigma'-1)))/((capn^alpha_1996) * (wages ^ (1-alpha_1996)))

*bys year sector: egen sum_tfpq_i_sigma = sum(tfpq_i^(`sigma'-1))
bys year sector: egen sum_tfpq_i_sigma = sum((tfpq_i/sum_count)^(`sigma'-1))
g a_sector_bar = sum_tfpq_i_sigma ^ (1/(`sigma'-1))

*g tfpq_i_deviation = tfpq_i * sum_count^(1/(`sigma'-1)) / a_sector_bar
g tfpq_i_deviation = tfpq_i  / a_sector_bar
g log_tfpq_i_deviation = log(tfpq_i_deviation)

 twoway (kdensity log_tfpq_i_deviation if year == 1996, ///
        xtitle("Distribution of TFPQ") ytitle("") title("Ghana: 1996 TFPQ of Manufacturing Firms")  ///
		note("Source: CSAE Ghana RPED/GMES Data" "Calculations as in Hsieh Klenow QJE 1996 paper") )
		
		graph export "$out/TFPQ Distribution 1996.pdf", as (pdf) replace


*Now recalculate TFPR Measures


g tfpr_i = routput / ( (wages^alpha_1996) * (capn^(1-alpha_1996)) )

*g tfpq_i= (routput^((`sigma')/(`sigma'-1)))/((capn^alpha_1996) * (wages ^ (1-alpha_1996)))


g tfpr_bar = sum_routput / ( (sum_wages^alpha_1996) * (sum_capn^(1-alpha_1996)) )

g log_tfpr_i = log(tfpr_i)
g log_tfpr_bar = log(tfpr_bar)

/* Alternate Measure: Calc is the Same 
g log_routput = log(routput)
g log_wages = log(wages)
g log_capn = log(capn)

g log_tfpr_i_alt = log_routput - (alpha_1996 * log_wages) - ((1-alpha_1996) * log_capn)
*/

g log_tfpr_i_deviation = log(tfpr_i/tfpr_bar)

 twoway (kdensity log_tfpr_i_deviation if year == 1996, ///
        xtitle("Distribution of TFPR") ytitle("") title("Ghana: 1996 TFPR of Manufacturing Firms")  ///
		note("Source: CSAE Ghana RPED/GMES Data" "Calculations as in Hsieh Klenow QJE 1996 paper") )
		
		graph export "$out/TFPR Distribution 1996.pdf", as (pdf) replace

		
		
*su log_tfpq, d

**Distribution

foreach prod in tfpr tfpq {
	bys year: egen sd_log_`prod' = sd(log_`prod'_i_deviation)
		foreach j in 10 25 75 90 {
			bys year: egen pct_`j'_log_`prod' = pctile(log_`prod'_i_deviation), p(`j')
		}
	g disp_75_25_log_`prod' = pct_75_log_`prod' - pct_25_log_`prod'
	g disp_90_10_log_`prod' = pct_90_log_`prod' - pct_10_log_`prod'
	
	bys year: egen nobs_`prod' = count(log_`prod'_i)
	}

preserve
g obs = [_n]
keep if obs == 1
keep year disp* sd_log* nobs_*
export excel "$out/Table I and II.xlsx", replace
restore
	

*Sources of Variation
encode(owndum), g(owndum_encode)
encode(locdum), g(locdum_encode)

label var owndum_encode "Ownership"
label var locdum_encode "Location"

save "$temp/All_Years_Log_TFPR_TFPQ", replace


foreach y in 1992 1996  1999 {
*foreach y in 1996 {
	use "$temp/All_Years_Log_TFPR_TFPQ", clear
	reg log_tfpr_i_deviation i.owndum_encode  if year == `y'
		regsave using "$temp/Sources of TFPR Variation in `y'", addlabel(Specification,"Ownership") replace
	reg log_tfpr_i_deviation i.owndum_encode fmage if year == `y'
		regsave using "$temp/Sources of TFPR Variation in `y'", addlabel(Specification,"+ Firm Age") append
	reg log_tfpr_i_deviation i.owndum_encode fmage worker skill workersq skillsq  if year == `y'
		regsave using "$temp/Sources of TFPR Variation in `y'", addlabel(Specification,"+ Workers") append
	reg log_tfpr_i_deviation i.owndum_encode fmage worker skill workersq skillsq i.locdum_encode if year == `y'
		regsave using "$temp/Sources of TFPR Variation in `y'", addlabel(Specification,"+ Region") append
	*reg log_tfpr_i_deviation i.owndum_encode fmage worker skill workersq skillsq eduwgt i.locdum_encode if year == `y'
	use "$temp/Sources of TFPR Variation in `y'", clear
	texsave using "$out/Underling Table III: Sources of TFPR Variation in `y'.tex", ///
	title(Regression for Sources of TFPR Variation) footnote("Source: CSAE Ghana RPED/GMES Data.  Percent sources represent R^2 from OLS regression adding each additional (set of) explanatory variable(s).") ///
	replace
	
	keep if var == "_cons"
	drop var
	foreach explan in Ownership Age Workers Region {
		g pre_`explan' = r2 if strpos(Specification,"`explan'") > 0
		egen `explan' = max(pre_`explan')
		drop pre_`explan'
		}
	g obs = [_n]
	keep if obs == 1
	keep Ownership Age Workers Region
	texsave using "$out/Table III: Sources of TFPR Variation in `y'.tex", ///
	title(Sources of TFPR Variation) footnote("Source: CSAE Ghana RPED/GMES Data.  Percent sources represent R^2 from OLS regression adding each additional (set of) explanatory variable(s).") ///
	replace
		
}
	

use "$temp/All_Years_Log_TFPR_TFPQ", clear

	 
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



