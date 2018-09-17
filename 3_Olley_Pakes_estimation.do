
set more off
clear

global raw "/Users/russellmorton/Desktop/SRA 2018/STATA/Data/Raw"

global temp "/Users/russellmorton/Desktop/SRA 2018/STATA/Data/Temp"

global cleaned "/Users/russellmorton/Desktop/SRA 2018/STATA/Data/Cleaned"

global out "/Users/russellmorton/Desktop/SRA 2018/STATA/Output/Olley Pakes Production"

use "$raw/firms", clear

*Check if there are any missing years


g drop = 0
replace drop = 1 if routput == . | routput == 0 | wages == . | wages == 0 | capn == . | capn == 0 | fmage == 0 

*drop if drop == 1

*Drop all observations 

*Naive Estimation: Unbalanced Panel

g lwages = ln(wages)

reg lroutput fmage lrcapn lwages
	regsave using "$temp/OlleyPakesResults", addlabel(Specification,"OLS") replace

xtreg lroutput fmage lrcapn lwages, fe
	regsave using "$temp/OlleyPakesResults", addlabel(Specification,"FE") append


*Naive Estiation: Balanced Panel

egen pre_min_year = min(year) if drop == 0
egen min_year = max(pre_min_year)
drop pre_min_year

egen pre_max_year = max(year) if drop == 0
egen max_year = max(pre_max_year)
drop pre_max_year

bys firm: egen pre_min_firm_year = min(year) if drop == 0
bys firm: egen min_firm_year = min(pre_min_firm_year) if drop == 0
drop pre_min_firm_year

bys firm: egen pre_max_firm_year = max(year) if drop == 0
bys firm: egen max_firm_year = max(pre_max_firm_year) if drop == 0
drop pre_max_firm_year

g balanced = min_firm_year == min_year & max_firm_year == max_year

reg lroutput fmage lrcapn lwages if balanced == 1 & drop == 0 
	regsave using "$temp/OlleyPakesResults", addlabel(Specification,"OLS_Balanced") append

xtreg lroutput fmage lrcapn lwages if balanced == 1 & drop == 0 , fe
	regsave using "$temp/OlleyPakesResults", addlabel(Specification,"FE_Balanced") append


///**** Olley Pakes Estimation ****///
*Stage 1

g linvest = log(rinv)
replace linvest = 0 if linvest == .

	forv i = 0(1)3 {
		forv j = 0(1)3 {
			forv k = 0(1)3 {
				capture g poly_linvest_`i'_lrcapn_`j'_fmage_`k' = linvest^`i' * lrcapn^`j'* fmage ^ `k'
			}
		}
	}

reg lroutput lwages poly* if drop == 0 

*Store estimates for the psi function and the coefficient on labor
g labor_beta_hat_stage1 = _b[lwages]
g capital_beta_hat_stage1 = _b[poly_linvest_0_lrcapn_1_fmage_0]
g fmage_beta_hat_stage1 = _b[poly_linvest_0_lrcapn_0_fmage_1]


g psi_hat_stage1 = 0
g beta_hat_use = 0

	forv i = 0(1)3 {
		forv j = 0(1)3 {
			forv k = 0(1)3 {
				replace beta_hat_use = _b[poly_linvest_`i'_lrcapn_`j'_fmage_`k']
				replace beta_hat_use = 0 if beta_hat_use == .
				replace psi_hat_stage1 = psi_hat_stage1 + beta_hat_use * poly_linvest_`i'_lrcapn_`j'_fmage_`k'
			}
		}
	}


*Stage 2: Selection Equation

xtset firm year
g survival = year < max_firm_year & drop == 0

*probit survival poly* i.firm i.year, iterate(50)
probit survival poly* i.year, iterate(500)
predict Psurvive_hat_probit

egen min_sector = min(sector)
egen max_sector = max(sector)
local minsector = min_firm
local maxsector = max_sector

forv s = `minsector'(1)`maxsector' {
	g sector_id_`s' = sector == `s'
}

local minyear = min_year
local maxyear = max_year

forv y = `minyear'(1)`maxyear' {
	g year_id_`y' = year == `y'
}

* gmm (survival - {xb: poly* firm_id* year_id*}), instruments(poly* firm_id* year_id*) i igmmit(500)
*gmm (survival - {xb: poly* year_id* sector_id*}), instruments(poly* year_id* sector_id*) i igmmit(500)
*predict Psurvive_hat_gmm
 
 xtset firm year
 
 g gfct_psi_hat = psi_hat_stage1 - L1.lrcapn * capital_beta_hat_stage1 - (fmage - 1) * fmage_beta_hat_stage1

	forv i = 0(1)3 {
		forv j = 0(1)3 {
			capture g poly_survive_probit_`i'_psi_`j' = Psurvive_hat_probit^`i' * gfct_psi_hat^`j'	
		}
	}
 
/*
foreach var in Psurvive_hat_gmm psi_hat_stage1 {
	forv i = 0(1)3 {
		forv j = 0(1)3 {
			capture g poly_survive_gmm_`i'_psi_`j' = Psurvive_hat_gmm^`i' * psi_hat_stage1^`j'	
		}
	}
}
*/
g output_less_labor = lroutput - labor_beta_hat_stage1 * lwages

save "$temp/OP dataset", replace

reg output_less_labor fmage lrcapn poly_survive_probit_1_psi_0 poly_survive_probit_2_psi_0 poly_survive_probit_3_psi_0 
		regsave using "$temp/OlleyPakesResults", addlabel(Specification,"OP_Only_P_OLS") append

xtreg output_less_labor fmage lrcapn poly_survive_probit_1_psi_0 poly_survive_probit_2_psi_0 poly_survive_probit_3_psi_0, fe
		regsave using "$temp/OlleyPakesResults", addlabel(Specification,"OP_Only_P_FE") append

reg output_less_labor fmage lrcapn poly_survive_probit_0_psi_1 poly_survive_probit_0_psi_2 poly_survive_probit_0_psi_3
		regsave using "$temp/OlleyPakesResults", addlabel(Specification,"OP_Only_Psi_OLS") append

xtreg output_less_labor fmage lrcapn poly_survive_probit_0_psi_1 poly_survive_probit_0_psi_2 poly_survive_probit_0_psi_3, fe
		regsave using "$temp/OlleyPakesResults", addlabel(Specification,"OP_Only_Psi_FE") append

reg output_less_labor fmage lrcapn poly_survive_probit_*
		regsave using "$temp/OlleyPakesResults", addlabel(Specification,"Semi_Parametric_OP_OLS") append

xtreg output_less_labor fmage lrcapn poly_survive_probit_*, fe
		regsave using "$temp/OlleyPakesResults", addlabel(Specification,"Semi_Parametric_OP_FE") append

		
		
*xtreg output_less_labor fmage lrcapn poly_survive_gmm_*, fe


***NEXT STEPS: 
*1.Make table
*2.Throw in lag of l and see if zero in last stage
*3.Productivity dispersion; compare to HK

use "$temp/OlleyPakesResults", clear

keep if inlist(var,"_cons","fmage","lrcapn","lwages") 

g obs = [_n]

g variable = "Labor" if obs == 1
replace variable = "Capital" if obs == 2
replace variable = "Age" if obs == 3
replace variable = "N" if obs == 4

levelsof Specification, local(specs)

g OP = strpos(Specification,"OP")

*di "specs is `specs'"

foreach sp in `specs' {
	di "sp is `sp'"


	g `sp' = 10000
	
	preserve
	keep if var == "lwages" & Specification == "`sp'"
	local value = coef
	local op = OP
	restore

	replace `sp' = `value' if variable == "Labor"

	preserve
	keep if var == "fmage" & Specification == "`sp'"
	local value = coef
	restore

	replace `sp' = `value' if variable == "Age"

	preserve
	keep if var == "lrcapn" & Specification == "`sp'"
	local value = coef
	restore

	replace `sp' = `value' if variable == "Capital"
	
	preserve
	keep if var == "lrcapn" & Specification == "`sp'"
	local value = N
	restore

	replace `sp' = `value' if variable == "N"
	
	replace `sp' = . if `op' >= 1 & variable == "Labor"

}

drop if obs > 4
drop var coef stderr N r2 Specification obs OP

texsave using "$out/OP Table IV: Estimates of Production Function Parameters.tex", ///
	title(Olley Pakes Production Function Regressions) footnote("Source: CSAE Ghana RPED/GMES Data") ///
	replace

	
/* LP Method */

use "$temp/OP dataset", clear

*Naive Estimates

xtset firm year

reg lroutput fmage lrcapn lwages lrmata
	regsave using "$temp/LPResults", addlabel(Specification,"OLS") replace

xtreg lroutput fmage lrcapn lwages, fe
	regsave using "$temp/LPResults", addlabel(Specification,"FE") append


reg lroutput fmage lrcapn lwages lrmata if balanced == 1 & drop == 0 
	regsave using "$temp/LPResults", addlabel(Specification,"OLS_Balanced") append

xtreg lroutput fmage lrcapn lwages lrmata if balanced == 1 & drop == 0 , fe
	regsave using "$temp/LPResults", addlabel(Specification,"FE_Balanced") append


*Stage 1

	forv i = 0(1)3 {
		forv j = 0(1)3 {
			forv k = 0(1)3 {
				capture g poly_lrmata_`i'_lrcapn_`j'_fmage_`k' = lrmata ^`i' * lrcapn^`j'* fmage ^ `k'
			}
		}
	}

reg lroutput lwages poly_lrmata* 

*Store estimates for the psi function and the coefficient on labor
g LP_labor_beta_hat_stage1 = _b[lwages]
g LP_capital_beta_hat_stage1 = _b[poly_lrmata_0_lrcapn_1_fmage_0]
g LP_fmage_beta_hat_stage1 = _b[poly_lrmata_0_lrcapn_0_fmage_1]


g psi_LP_hat_stage1 = 0
g beta_LP_hat_use = 0

	forv i = 0(1)3 {
		forv j = 0(1)3 {
			forv k = 0(1)3 {
				replace beta_LP_hat_use = 0
				capture replace beta_LP_hat_use = _b[poly_lrmata_`i'_lrcapn_`j'_fmages_`k']
				replace beta_LP_hat_use = 0 if beta_hat_use == .
				replace psi_LP_hat_stage1 = psi_LP_hat_stage1 + beta_LP_hat_use * poly_lrmata_`i'_lrcapn_`j'_fmage_`k'
			}
		}
	}

	
probit survival poly_lrmata* i.year, iterate(500)
predict Psurvive_LP_hat_probit

 xtset firm year
 
 g LP_gfct_psi_hat = psi_LP_hat_stage1 - L1.lrcapn * LP_capital_beta_hat_stage1 - (fmage - 1) * LP_fmage_beta_hat_stage1


	forv i = 0(1)3 {
		forv j = 0(1)3 {
			capture g LP_poly_survive_probit_`i'_psi_`j' = Psurvive_LP_hat_probit^`i' * LP_gfct_psi_hat ^`j'	
		}
	}
	
g output_LP_less_labor = lroutput - LP_labor_beta_hat_stage1 * lwages


save "$temp/LP dataset", replace

reg output_LP_less_labor fmage lrcapn lrmata LP_poly_survive_probit_1_psi_0 LP_poly_survive_probit_2_psi_0 LP_poly_survive_probit_3_psi_0 
		regsave using "$temp/LPResults", addlabel(Specification,"LP_Only_P_OLS") append

xtreg output_LP_less_labor fmage lrcapn lrmata LP_poly_survive_probit_1_psi_0 LP_poly_survive_probit_2_psi_0 LP_poly_survive_probit_3_psi_0, fe
		regsave using "$temp/LPResults", addlabel(Specification,"LP_Only_P_FE") append

reg output_LP_less_labor fmage lrcapn lrmata LP_poly_survive_probit_0_psi_1 LP_poly_survive_probit_0_psi_2 LP_poly_survive_probit_0_psi_3
		regsave using "$temp/LPResults", addlabel(Specification,"LP_Only_Psi_OLS") append

xtreg output_LP_less_labor fmage lrcapn lrmata LP_poly_survive_probit_0_psi_1 LP_poly_survive_probit_0_psi_2 LP_poly_survive_probit_0_psi_3, fe
		regsave using "$temp/LPResults", addlabel(Specification,"LP_Only_Psi_FE") append

reg output_LP_less_labor fmage lrcapn lrmata LP_poly_survive_probit_*
		regsave using "$temp/LPResults", addlabel(Specification,"Semi_Parametric_LP_OLS") append

xtreg output_LP_less_labor fmage lrcapn lrmata LP_poly_survive_probit_*, fe
		regsave using "$temp/LPResults", addlabel(Specification,"Semi_Parametric_LP_FE") append


		
use "$temp/LPResults", clear

keep if inlist(var,"_cons","fmage","lrcapn","lwages","lrmata") 

g obs = [_n]

g variable = "Labor" if obs == 1
replace variable = "Capital" if obs == 2
replace variable = "Age" if obs == 3
replace variable = "Materials" if obs == 4
replace variable = "N" if obs == 5

levelsof Specification, local(specs)

g LP = strpos(Specification,"LP")

*di "specs is `specs'"

foreach sp in `specs' {
	di "sp is `sp'"

	g `sp' = 10000
	
	preserve
	keep if var == "lwages" & Specification == "`sp'"
	local value = coef
	local lp = LP
	restore

	replace `sp' = `value' if variable == "Labor"

	preserve
	keep if var == "fmage" & Specification == "`sp'"
	local value = coef
	restore

	replace `sp' = `value' if variable == "Age"

	preserve
	keep if var == "lrcapn" & Specification == "`sp'"
	local value = coef
	restore

	replace `sp' = `value' if variable == "Capital"
	
	preserve
	keep if var == "lrcapn" & Specification == "`sp'"
	local value = N
	restore

	replace `sp' = `value' if variable == "N"
	
	preserve
	keep if var == "lrmata" & Specification == "`sp'"
	local value = coef
	restore

	replace `sp' = `value' if variable == "Materials"
	
	
	replace `sp' = . if `lp' > 0 & variable == "Labor"

}

drop if obs > 5
drop var coef stderr N r2 Specification obs LP

*Reorder Columns


order variable FE FE_Balanced OLS OLS_Balanced LP_Only_P_FE LP_Only_P_OLS LP_Only_Psi_FE LP_Only_Psi_OLS Semi_Parametric_LP_FE Semi_Parametric_LP_OLS


texsave using "$out/LP Table IV: Estimates of Production Function Parameters.tex", ///
	title(Levinsohn Petrin Production Function Regressions) footnote("Source: CSAE Ghana RPED/GMES Data") ///
	replace

/*******
Ackerberg Method
*******/

use "$temp/OP dataset", clear

	forv i = 0(1)3 {
		forv j = 0(1)3 {
			forv k = 0(1)3 {
				forv l = 0(1)3 {
					g py_mat_`i'_cap_`j'_fmage_`k'_labor_`l' = lrmata ^`i' * lrcapn^`j'* fmage ^ `k' * lwages ^ `l'
				}
			}
		}
	}
 
reg lroutput py_mat* 

g AC_labor_beta_hat_stage1 = _b[py_mat_0_cap_0_fmage_0_labor_1]
g AC_capital_beta_hat_stage1 = _b[py_mat_0_cap_1_fmage_0_labor_0]
g AC_fmage_beta_hat_stage1 = _b[py_mat_0_cap_0_fmage_1_labor_0]


g psi_AC_hat_stage1 = 0
g beta_AC_hat_use = 0

	forv i = 0(1)3 {
		forv j = 0(1)3 {
			forv k = 0(1)3 {
				forv l = 0(1)3 {
					replace beta_AC_hat_use = 0
					capture replace beta_AC_hat_use = _b[py_mat_`i'_cap_`j'_fmage_`k'_labor_`l']
					replace beta_AC_hat_use = 0 if beta_AC_hat_use == .
					replace psi_AC_hat_stage1 = psi_AC_hat_stage1 + beta_AC_hat_use * py_mat_`i'_cap_`j'_fmage_`k'_labor_`l'
			}
		}
	}
}

	
probit survival py_mat* i.year, iterate(500)
predict Psurvive_AC_hat_probit

 xtset firm year
 
 g AC_gfct_psi_hat = psi_AC_hat_stage1 - L1.lrcapn * AC_capital_beta_hat_stage1 - (fmage - 1) * AC_fmage_beta_hat_stage1 ///
							- L1.lwages * AC_labor_beta_hat_stage1 


	forv i = 0(1)3 {
		forv j = 0(1)3 {
			capture g AC_poly_survive_probit_`i'_psi_`j' = Psurvive_AC_hat_probit^`i' * AC_gfct_psi_hat ^`j'	
		}
	}
	

save "$temp/AC dataset", replace

reg routput fmage lrcapn lrmata lwages AC_poly_survive_probit_1_psi_0 AC_poly_survive_probit_2_psi_0 AC_poly_survive_probit_3_psi_0 
		regsave using "$temp/ACResults", addlabel(Specification,"AC_Only_P_OLS") replace

xtreg routput fmage lrcapn lrmata lwages AC_poly_survive_probit_1_psi_0 AC_poly_survive_probit_2_psi_0 AC_poly_survive_probit_3_psi_0, fe
		regsave using "$temp/ACResults", addlabel(Specification,"AC_Only_P_FE") append

reg routput fmage lrcapn lrmata lwages AC_poly_survive_probit_0_psi_1 AC_poly_survive_probit_0_psi_2 AC_poly_survive_probit_0_psi_3
		regsave using "$temp/ACResults", addlabel(Specification,"AC_Only_Psi_OLS") append

xtreg routput fmage lrcapn lrmata lwages AC_poly_survive_probit_0_psi_1 AC_poly_survive_probit_0_psi_2 AC_poly_survive_probit_0_psi_3, fe
		regsave using "$temp/ACResults", addlabel(Specification,"AC_Only_Psi_FE") append

reg routput fmage lrcapn lrmata lwages AC_poly_survive_probit_*
		regsave using "$temp/ACResults", addlabel(Specification,"Semi_Parametric_AC_OLS") append

xtreg routput fmage lrcapn lrmata lwages AC_poly_survive_probit_*, fe
		regsave using "$temp/ACResults", addlabel(Specification,"Semi_Parametric_AC_FE") append

		
use "$temp/ACResults", clear

keep if inlist(var,"_cons","fmage","lrcapn","lwages","lrmata") 

g obs = [_n]

g variable = "Labor" if obs == 1
replace variable = "Capital" if obs == 2
replace variable = "Age" if obs == 3
replace variable = "Materials" if obs == 4
replace variable = "N" if obs == 5

levelsof Specification, local(specs)

g LP = strpos(Specification,"LP")

*di "specs is `specs'"

foreach sp in `specs' {
	di "sp is `sp'"

	g `sp' = 10000
	
	preserve
	keep if var == "lwages" & Specification == "`sp'"
	local value = coef
	restore

	replace `sp' = `value' if variable == "Labor"

	preserve
	keep if var == "fmage" & Specification == "`sp'"
	local value = coef
	restore

	replace `sp' = `value' if variable == "Age"

	preserve
	keep if var == "lrcapn" & Specification == "`sp'"
	local value = coef
	restore

	replace `sp' = `value' if variable == "Capital"
	
	preserve
	keep if var == "lrcapn" & Specification == "`sp'"
	local value = N
	restore

	replace `sp' = `value' if variable == "N"
	
	preserve
	keep if var == "lrmata" & Specification == "`sp'"
	local value = coef
	restore

	replace `sp' = `value' if variable == "Materials"
	
	
}

drop if obs > 5

drop var coef stderr N r2
drop Specification obs
drop LP


*order variable FE FE_Balanced OLS OLS_Balanced AC_Only_P_FE AC_Only_P_OLS AC_Only_Psi_FE AC_Only_Psi_OLS Semi_Parametric_AC_FE Semi_Parametric_AC_OLS
order variable AC_Only_P_FE AC_Only_P_OLS AC_Only_Psi_FE AC_Only_Psi_OLS Semi_Parametric_AC_FE Semi_Parametric_AC_OLS


texsave using "$out/Ack Table IV: Estimates of Production Function Parameters.tex", ///
	title(Ackerberg Production Function Regressions) footnote("Source: CSAE Ghana RPED/GMES Data") ///
	replace
