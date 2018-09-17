*Created by RM on 2018.05.20 for Dynamic Inputs and Resource Mis(Allocation) style Works

set more off
clear

global raw "/Users/russellmorton/Desktop/SRA 2018/STATA/Data/Raw"

global temp "/Users/russellmorton/Desktop/SRA 2018/STATA/Data/Temp"

global cleaned "/Users/russellmorton/Desktop/SRA 2018/STATA/Data/Cleaned"

global out "/Users/russellmorton/Desktop/SRA 2018/STATA/Output/Dynamic Inputs Misallocation"

use "$raw/firms", clear

///***
*Begin by estimating betas for each industry
***///

 drop if routput == . | routput == 0 | wages == . | wages == 0 | capn == . | capn == 0 | fmage == 0 | capn < 0

 *As in previous samples use 1996
 
 g labor_share_1996 = wages / routput if year == 1996
 g materials_share_1996 = rmata / routput if year == 1996
 
 bys sector: egen pre_median_labor_share_1996 = median(labor_share_1996) if labor_share_1996 != .
 bys sector: egen pre_median_materials_share_1996 = median(materials_share_1996) if materials_share_1996 != .
 
 bys sector: egen median_labor_share_1996 =  max(pre_median_labor_share_1996)
 bys sector: egen median_materials_share_1996 =  max(pre_median_materials_share_1996)
 
 g beta_labor_1996 = median_labor_share_1996
 g beta_materials_1996 = median_materials_share_1996
 g beta_cap_1996 = (4/3) - beta_labor_1996 - beta_materials_1996
 
 
 g lwages = ln(wages)
 g tfpr = lroutput - beta_cap_1996 * lcap - beta_labor_1996 * lwages - beta_materials_1996 * lrmata
 
 drop if tfpr == .
 
 g mrpk = ln(beta_cap_1996) + lroutput - lcap
 
* g linvest = log(rinv)

///***
*Reduced Form Results
***///

**Table 2:
*Medians

egen median_workers = median(worker)

xtset firm year

egen median_delta_s = median(D1.lroutput)
egen median_delta_tfpr = median(D1.tfpr)

*Standard Deviations

egen sd_mrpk = sd(mrpk)
egen sd_k = sd(lcap)
egen sd_tfpr = sd(tfpr)
xtset firm year
g change_tfpr = D1.tfpr
egen vol = sd(change_tfpr)

**Figure 2
bys sector: egen ind_sd_mrpk = sd(mrpk)
bys sector: egen ind_vol = sd(change_tfpr)

egen tag_sector = tag(sector)

preserve
keep if tag_sector == 1

twoway (scatter ind_sd_mrpk  ind_vol) (lfit ind_sd_mrpk ind_vol)
reg ind_sd_mrpk ind_vol
reg ind_sd_mrpk ind_vol i.year

restore

**Figure 2 Alt: Use 1996 Alone

bys sector: egen pre_ind_sd_mrpk_1996 = sd(mrpk) if year == 1996
bys sector: egen ind_sd_mrpk_1996 = max(pre_ind_sd_mrpk_1996)
bys sector: egen pre_ind_vol_1996 = sd(change_tfpr) if year == 1996
bys sector: egen ind_vol_1996 = max(pre_ind_vol_1996)

preserve
keep if year == 1996
egen tag_sector_1996 = tag(sector)
keep if tag_sector_1996 == 1

twoway (scatter ind_sd_mrpk_1996  ind_vol_1996) (lfit ind_sd_mrpk_1996 ind_vol_1996)
reg ind_sd_mrpk_1996 ind_vol_1996
restore 

**Table 3
bys sector year: egen ind_year_sd_mrpk = sd(mrpk)
bys sector year: egen ind_year_vol = sd(change_tfpr)

egen tag_sector_year = tag(sector year)

reg ind_year_sd_mrpk ind_year_vol i.year if tag_sector_year == 1, cluster(sector) 


**Table 4: Robustness

*col 1:
reg ind_year_sd_mrpk ind_vol if tag_sector_year == 1, cluster(sector)
*col 1 alt:
reg ind_year_sd_mrpk ind_vol i.year if tag_sector_year == 1, cluster(sector)

*col 2:
xtset firm year
g lag_tfpr = L1.tfpr

egen minsector = min(sector)
egen maxsector = max(sector)

local minsec = minsector
local maxsec = maxsector

g ar_tfpr_rho_hat = .
g ar_tfpr_mu_hat = .
*g tfpr_ar_resid = .

forv i = `minsec'(1)`maxsec' {
	reg tfpr lag_tfpr if sector == `i'
	replace ar_tfpr_rho_hat = _b[lag_tfpr] if sector == `i'
	replace ar_tfpr_mu_hat =  _b[_cons] if sector == `i'	
	*predict tfpr_ar_resid_pred, resid 
	*replace tfpr_ar_resid = tfpr_ar_resid_pred if sector == `i'
	*drop tfpr_ar_resid_pred
}

g tfpr_ar_residual = tfpr - ar_tfpr_mu_hat - (ar_tfpr_rho_hat * lag_tfpr)
g tfpr_ar_residual_2 = tfpr_ar_residual^2

bys sector: egen count_sector = count(sector) if tfpr_ar_residual != .
bys sector: egen sum_tfpr_ar_resid_2 = sum(tfpr_ar_residual_2)

g ar_tfpr_sigma_hat = sum_tfpr_ar_resid_2 / count_sector

reg ind_year_sd_mrpk ar_tfpr_sigma_hat

*col 3:

g ar_tfpr_fe_rho_hat = .
g ar_tfpr_fe_mu_hat = .
*g tfpr_ar_resid = .

forv i = `minsec'(1)`maxsec' {
	reg tfpr lag_tfpr i.firm if sector == `i'
	replace ar_tfpr_fe_rho_hat = _b[lag_tfpr] if sector == `i'
	replace ar_tfpr_fe_mu_hat =  _b[_cons] if sector == `i'	
	*predict tfpr_ar_resid_pred, resid 
	*replace tfpr_ar_resid = tfpr_ar_resid_pred if sector == `i'
	*drop tfpr_ar_resid_pred
}

g tfpr_ar_fe_residual = tfpr - ar_tfpr_fe_mu_hat - (ar_tfpr_fe_rho_hat * lag_tfpr)
g tfpr_ar_fe_residual_2 = tfpr_ar_fe_residual^2

bys sector: egen count_sector_fe = count(sector) if tfpr_ar_fe_residual != .
bys sector: egen sum_tfpr_ar_fe_resid_2 = sum(tfpr_ar_fe_residual_2)

g ar_tfpr_fe_sigma_hat = sum_tfpr_ar_fe_resid_2 / count_sector_fe

reg ind_year_sd_mrpk ar_tfpr_fe_sigma_hat


*Table 5
xtset firm year
reg mrpk D1.tfpr lcap L1.tfpr i.year i.sector


*Table 6

bys sector year: egen median_vol_ind_year = median(ind_year_vol)
g vol_over_median_dummy = ind_year_vol > median_vol_ind_year
g vol_over_median_interaction = ind_year_vol * vol_over_median_dummy

*section 1:

xtset firm year
g delta_mrpk = D1.mrpk
bys year sector: egen sd_delta_mrpk = sd(delta_mrpk)
reg sd_delta_mrpk ind_year_vol if tag_sector_year == 1

*section 2:
xtset firm year
g delta_cap = D1.lcap
bys year sector: egen sd_delta_cap = sd(delta_cap)
reg sd_delta_cap ind_year_vol if tag_sector_year == 1

*section 3:
reg sd_delta_cap ind_year_vol vol_over_median_interaction if tag_sector_year == 1


*Table 7

*col 1:
bys year sector: egen sd_year_ind_mrpk = sd(mrpk)
*col 2:
g mrpl = ln(beta_labor_1996) + lroutput - lwages
bys year sector: egen sd_year_ind_mrpl = sd(mrpl)
*col 3:
g mrpm = ln(beta_materials_1996) + lroutput - lrmata
bys year sector: egen sd_year_ind_mrpm = sd(mrpm)


/*
*STRUCTURAL ANALYSIS 
*/

*Estimate the adjustment cost function
*Moments:
*1. % of firms with less than 5% year-on-year change in capital
xtset firm year
g change_cap_populated = 0
replace change_cap_populated = 1 if D1.cap != .
g change_cap_perc = abs(D1.cap) / L1.cap
g change_cap_perc_lt_5 = 0 if change_cap_populated == 1
replace change_cap_perc_lt_5 = 1 if change_cap_perc < .05

bys sector: egen count_cap_change_pop = sum(change_cap_populated)
bys sector: egen count_cap_change_lt_5_perc = sum(change_cap_perc_lt_5)
g moment_cap_change_lt_5_perc = count_cap_change_lt_5_perc / count_cap_change_pop

*2. % of firms with more than 20% year-on-year change in capital
g change_cap_perc_gt_20 = 0 if change_cap_populated == 1
replace change_cap_perc_gt_20 = 1 if change_cap_perc > .2

bys sector: egen count_cap_change_gt_20_perc = sum(change_cap_perc_gt_20)
g moment_cap_change_gt_20_perc = count_cap_change_gt_20_perc / count_cap_change_pop

*3. SD of change in log capital
xtset firm year
g log_cap_change = D1.lcap
bys sector: egen moment_sd_change_cap = sd(log_cap_change)

save "$temp/StructuralAnalysisData", replace

**Need to compute optimal policy for firms as function of the fixed and variable cost given omega and capital




/*

*MLE estimator
g sigma1_best = -10
g rho1_best = -1000
g mu1_best = -100000000
g sum_likelihood1_best = -10000000000000

*local i = 0

keep if sector == 3

forv sigma1 = 1(10)501 {
	*forv rho1 = -20(2)20 {
		forv rho1 = -5(.25)5 {
			forv mu1 = -300(10)300 {
		
			*local i = `i' + 1
			
			xtset firm year
			g likelihood1 = 0
			replace likelihood1 = -(tfpr-`rho1'*L1.tfpr - `mu1')^2 / (2*`sigma1'^2)
			egen sum_likelihood1 = sum(likelihood1)	
			replace sigma1_best = `sigma1' if sum_likelihood1 > sum_likelihood1_best
			replace rho1_best = `rho1' if sum_likelihood1 > sum_likelihood1_best
			replace mu1_best = `mu1' if sum_likelihood1 > sum_likelihood1_best
			replace sum_likelihood1_best = sum_likelihood1 if sum_likelihood1 > sum_likelihood1_best
			
			su sum_likelihood1
			su sum_likelihood1_best
			
			*g sum_likelihood_sigma_`sigma1'_rho_`rho1'_mu_`mu1' = sum_likelihood1

			drop likelihood1 sum_likelihood1
			
			*g likelihood_best_`i' = sum_likelihood1_best
		}
	}
}

	
su sigma1_best
su rho1_best	
su mu1_best

*/





