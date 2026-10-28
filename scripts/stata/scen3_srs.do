version 19
set more off
set linesize 255

* Scenario 3 without the survey design: a simple random sample of households drawn from the population,
* all members kept, fitted unweighted. Separates the data generation from the weighting and design.
* Usage: do scen3_srs.do <sim_id> <population>
*   population: dec2024 or oct2026
* The sample is self-weighting, so the one unweighted fit replaces both scen3 steps.
* SEs are clustered on households.
args sim_id population

if !inlist("`population'", "dec2024", "oct2026") {
    display as error "population must be dec2024 or oct2026"
    exit 198
}
local repo "/work/users/a/q/aquijano/mixed_sims/hchs_mixedmodels_simulations"
local population_file "/work/users/b/e/beibo/v3data/Dec2024/population_3visits_Dec2024.csv"
if "`population'" == "oct2026" {
    local population_file "/work/users/b/e/beibo/v3data/Oct2026/population_3visits_Oct2026.csv"
}
local out_dir "`repo'/data/estimates/`population'/scen3_srs"
capture mkdir "`repo'/data/estimates"
capture mkdir "`repo'/data/estimates/`population'"
capture mkdir "`out_dir'"

* About the number of households in a design sample
local n_hh 11500

import delimited "`population_file'", clear varnames(1)
rename *, lower
keep bgid hhid subid age_base x6 x12 x17 x18 y_bmi y_gfr

egen long gl_hhid  = group(bgid hhid)
egen long gl_subid = group(bgid hhid subid)

preserve
keep gl_hhid
duplicates drop
sort gl_hhid
set seed `=200000 + `sim_id''
gen double u = runiform()
sort u
keep in 1/`n_hh'
keep gl_hhid
tempfile chosen
save "`chosen'"
restore
merge m:1 gl_hhid using "`chosen'", keep(match) nogenerate

gen byte age_strat_new = (age_base >= 45)

meglm y_gfr age_strat_new x12 x17 x18 x6 y_bmi || gl_hhid: || gl_subid: , vce(cluster gl_hhid)
display "FITCHECK sim=`sim_id' converged=" e(converged)

tempname memhold rt
tempfile results
postfile `memhold' str32 term double estimate double std_error double conf_low double conf_high ///
    using "`results'", replace

matrix `rt' = r(table)
foreach var in _cons age_strat_new x12 x17 x18 x6 y_bmi {
    local term = cond("`var'" == "_cons", "(Intercept)", "`var'")
    local col = colnumb(`rt', "y_gfr:`var'")
    post `memhold' ("`term'") (`rt'[1, `col']) (`rt'[2, `col']) (`rt'[5, `col']) (`rt'[6, `col'])
}
foreach pair in "var_hhid /:var(_cons[gl_hhid])" "var_subid /:var(_cons[gl_hhid>gl_subid])" "s2 /:var(e.y_gfr)" {
    gettoken term name : pair
    local col = colnumb(`rt', trim("`name'"))
    post `memhold' ("`term'") (`rt'[1, `col']) (`rt'[2, `col']) (`rt'[5, `col']) (`rt'[6, `col'])
}
postclose `memhold'

use "`results'", clear
export delimited using "`out_dir'/scen3_srs_coef_`sim_id'.csv", replace
