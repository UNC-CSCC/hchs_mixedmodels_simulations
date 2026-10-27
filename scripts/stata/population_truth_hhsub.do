version 19
set more off
set linesize 255

* Census ("true") parameters: the analysis model fitted, unweighted, to the whole population with
* household and subject random intercepts. Gives the fixed effects plus var_hhid, var_subid and s2
* used as the truth for coverage of scen9/scen15/scen16.
* Usage: do population_truth_hhsub.do <population> [last_row]
*   population: dec2024 or oct2026
*   last_row:   optional; fit only the first <last_row> rows (smoke test). Omit for the full population.
args population last_row

if !inlist("`population'", "dec2024", "oct2026") {
    display as error "population must be dec2024 or oct2026"
    exit 198
}
local repo "/work/users/a/q/aquijano/mixed_sims/hchs_mixedmodels_simulations"
local population_file "/work/users/b/e/beibo/v3data/Dec2024/population_3visits_Dec2024.csv"
if "`population'" == "oct2026" {
    local population_file "/work/users/b/e/beibo/v3data/Oct2026/population_3visits_Oct2026.csv"
}
capture mkdir "`repo'/data/population_coefficients"
local out_file "`repo'/data/population_coefficients/population_coef_`population'_hhsub.csv"
if "`last_row'" != "" {
    local row_range "rowrange(1:`last_row')"
    local out_file "`repo'/data/population_coefficients/smoke_population_coef_`population'_hhsub.csv"
}

import delimited "`population_file'", clear varnames(1) `row_range'
rename *, lower

gen byte age_strat_new = (age_base >= 45)
egen long gl_hhid  = group(bgid hhid)
egen long gl_subid = group(bgid hhid subid)

mixed y_gfr age_strat_new x12 x17 x18 x6 y_bmi || gl_hhid: || gl_subid: , mle
display "FITCHECK population converged=" e(converged)

tempname memhold
tempfile results
postfile `memhold' str32 term double estimate using "`results'", replace
foreach var in _cons age_strat_new x12 x17 x18 x6 y_bmi {
    local term = cond("`var'" == "_cons", "(Intercept)", "`var'")
    post `memhold' ("`term'") (_b[`var'])
}
* mixed estimates log standard deviations: lns1 = household level, lns2 = subject level, lnsig_e = residual
post `memhold' ("var_hhid")  (exp(2 * _b[lns1_1_1:_cons]))
post `memhold' ("var_subid") (exp(2 * _b[lns2_1_1:_cons]))
post `memhold' ("s2")        (exp(2 * _b[lnsig_e:_cons]))
postclose `memhold'

use "`results'", clear
export delimited using "`out_file'", replace
