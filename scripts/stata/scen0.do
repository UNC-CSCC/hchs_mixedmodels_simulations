version 19
set more off
set linesize 255

* Scenario 0: one-step pseudo-likelihood for the 3-level model (visits in subjects in households)
* Usage: do scen0.do <sim_id> <population>
*   population: dec2024 (v3data/Dec2024 samples) or oct2026 (v3data/Oct2026 samples)
* Household (upper level): unconditional weight w_bghh, unscaled.
* Subject: w_sub, scaled globally to mean 1 across the sample (no within-household scaling).
* SEs clustered on households.
args sim_id population

if !inlist("`population'", "dec2024", "oct2026") {
    display as error "population must be dec2024 or oct2026"
    exit 198
}
local repo "/work/users/a/q/aquijano/mixed_sims/hchs_mixedmodels_simulations"
local data_dir "/work/users/b/e/beibo/v3data/Dec2024/sample_miss00pct"
if "`population'" == "oct2026" {
    local data_dir "/work/users/b/e/beibo/v3data/Oct2026/sample_miss00pct"
}
local out_dir "`repo'/data/estimates/`population'/scen0"
capture mkdir "`repo'/data/estimates"
capture mkdir "`repo'/data/estimates/`population'"
capture mkdir "`out_dir'"

import delimited "`data_dir'/samplemiss_`sim_id'.csv", clear varnames(1) bindquote(strict)
rename *, lower

gen byte age_strat_new = (age_base >= 45)
egen long gl_hhid  = group(bgid hhid)
egen long gl_subid = group(bgid hhid subid)

keep if !missing(y_gfr, x12, x6, x17, y_bmi, x18, x14, age_strat_new, bgid, hhid, subid, w_bghh, w_sub) ///
      & w_bghh > 0 & w_sub > 0

gen w_visit = 1
egen double sum_w_sub = total(w_sub)
gen double w_sub_glob = w_sub * (_N / sum_w_sub)

meglm y_gfr age_strat_new x12 x17 x18 x6 y_bmi [pw=w_visit] || ///
    gl_hhid:  , pweight(w_bghh) || ///
    gl_subid: , pweight(w_sub_glob) ///
    vce(cluster gl_hhid)
display "FITCHECK sim=`sim_id' converged=" e(converged)

tempname memhold rt
tempfile results
postfile `memhold' str32 term double estimate double std_error double conf_low double conf_high ///
    using "`results'", replace

matrix `rt' = r(table)
foreach var in _cons age_strat_new x12 x17 x18 x6 y_bmi {
    local col  = colnumb(`rt', "y_gfr:`var'")
    local term = cond("`var'" == "_cons", "(Intercept)", "`var'")
    post `memhold' ("`term'") (`rt'[1, `col']) (`rt'[2, `col']) (`rt'[5, `col']) (`rt'[6, `col'])
}
foreach pair in "var_hhid /:var(_cons[gl_hhid])" "var_subid /:var(_cons[gl_hhid>gl_subid])" "s2 /:var(e.y_gfr)" {
    gettoken term name : pair
    local col = colnumb(`rt', trim("`name'"))
    post `memhold' ("`term'") (`rt'[1, `col']) (`rt'[2, `col']) (`rt'[5, `col']) (`rt'[6, `col'])
}
postclose `memhold'

use "`results'", clear
export delimited using "`out_dir'/scen0_coef_`sim_id'.csv", replace
