version 19
set more off
set linesize 255

* Scenario 3 model fitted to the design samples without weights. Compared with scen3 (weighted) it shows
* whether the design is informative for the variance components or the weighting adds the bias.
* Usage: do scen3_unw.do <sim_id> <population>
*   population: dec2024 (v3data/Dec2024 samples) or oct2026 (v3data/Oct2026 samples)
* Without weights both scen3 steps are the same fit, so it is fitted once. SEs are clustered on households.
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
local out_dir "`repo'/data/estimates/`population'/scen3_unw"
capture mkdir "`repo'/data/estimates"
capture mkdir "`repo'/data/estimates/`population'"
capture mkdir "`out_dir'"

import delimited "`data_dir'/samplemiss_`sim_id'.csv", clear varnames(1) bindquote(strict)
rename *, lower

gen byte age_strat_new = (age_base >= 45)
egen long gl_hhid  = group(bgid hhid)
egen long gl_subid = group(bgid hhid subid)

* Same records as scen3
keep if !missing(y_gfr, x12, x6, x17, y_bmi, x18, x14, age_strat_new, bgid, hhid, subid, w_sub, w_bghh, sub_s1) ///
      & w_sub > 0 & w_bghh > 0 & sub_s1 > 0

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
export delimited using "`out_dir'/scen3_unw_coef_`sim_id'.csv", replace
