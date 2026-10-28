version 19
set more off
set linesize 255

* Scenario 3 with design-stratum indicators (strat 2-8) added to the model in both steps.
* Usage: do scen3_strata.do <sim_id> <population>
*   population: dec2024 (v3data/Dec2024 samples) or oct2026 (v3data/Oct2026 samples)
* Step 1 (globally scaled weights): fixed effects and residual variance.
* Step 2 (unscaled household weight w_bghh, effective-size subject weight sub_s1): household and subject
* variances, with all fixed effects (strata included) and the residual variance held at Step 1 values.
* Both steps cluster SEs on households.
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
local out_dir "`repo'/data/estimates/`population'/scen3_strata"
capture mkdir "`repo'/data/estimates"
capture mkdir "`repo'/data/estimates/`population'"
capture mkdir "`out_dir'"

import delimited "`data_dir'/samplemiss_`sim_id'.csv", clear varnames(1) bindquote(strict)
rename *, lower

gen byte age_strat_new = (age_base >= 45)
egen long gl_hhid  = group(bgid hhid)
egen long gl_subid = group(bgid hhid subid)

local strata ""
forvalues k = 2/8 {
    gen byte strat_`k' = (strat == `k')
    local strata "`strata' strat_`k'"
}
local covars "age_strat_new x12 x17 x18 x6 y_bmi `strata'"

keep if !missing(y_gfr, x12, x6, x17, y_bmi, x18, x14, age_strat_new, strat, bgid, hhid, subid, w_sub, w_bghh, sub_s1) ///
      & w_sub > 0 & w_bghh > 0 & sub_s1 > 0

gen w_visit = 1

* Step 1: globally scaled subject weights (mean 1); fixed effects and s2
egen double sum_w_sub_total  = total(w_sub)
gen  double w_sub_glob       = w_sub * (_N / sum_w_sub_total)

meglm y_gfr `covars' [pw=w_visit] || ///
    gl_hhid:  , pweight(w_bghh) || ///
    gl_subid: , pweight(w_sub_glob)  ///
    vce(cluster gl_hhid)
display "FITCHECK sim=`sim_id' step=1 converged=" e(converged)

tempname rt1
matrix `rt1' = r(table)

foreach var in _cons `covars' {
    local b_`var'  = _b[`var']
    local se_`var' = _se[`var']
}

local col_s2_step1 = colnumb(`rt1', "/:var(e.y_gfr)")
local s2_step1     = `rt1'[1, `col_s2_step1']
local se_s2_step1  = `rt1'[2, `col_s2_step1']
local ll_s2_step1  = `rt1'[5, `col_s2_step1']
local ul_s2_step1  = `rt1'[6, `col_s2_step1']

* Step 2: fixed effects and s2 held at Step 1; household and subject variances
constraint drop _all
local c = 0
foreach var in _cons `covars' {
    local c = `c' + 1
    constraint `c' [y_gfr]`var' = `b_`var''
}
local c = `c' + 1
constraint `c' [/]var(e.y_gfr) = `s2_step1'

meglm y_gfr `covars' [pw=w_visit] || ///
    gl_hhid:  , pweight(w_bghh) || ///
    gl_subid: , pweight(sub_s1) ///
    vce(cluster gl_hhid) ///
    constraints(1/`c')
display "FITCHECK sim=`sim_id' step=2 converged=" e(converged)

tempname memhold rt2
tempfile results
postfile `memhold' str32 term double estimate double std_error double conf_low double conf_high ///
    using "`results'", replace

matrix `rt2' = r(table)

local col_hh  = colnumb(`rt2', "/:var(_cons[gl_hhid])")
local var_hh  = `rt2'[1, `col_hh']
local se_hh   = `rt2'[2, `col_hh']
local ll_hh   = `rt2'[5, `col_hh']
local ul_hh   = `rt2'[6, `col_hh']

local col_sub = colnumb(`rt2', "/:var(_cons[gl_hhid>gl_subid])")
local var_sub = `rt2'[1, `col_sub']
local se_sub  = `rt2'[2, `col_sub']
local ll_sub  = `rt2'[5, `col_sub']
local ul_sub  = `rt2'[6, `col_sub']

local z975 = invnormal(0.975)
foreach var in _cons `covars' {
    local term = cond("`var'" == "_cons", "(Intercept)", "`var'")
    post `memhold' ("`term'") (`b_`var'') (`se_`var'') ///
        (`b_`var'' - `z975'*`se_`var'') (`b_`var'' + `z975'*`se_`var'')
}

post `memhold' ("var_hhid")  (`var_hh')  (`se_hh')  (`ll_hh')  (`ul_hh')
post `memhold' ("var_subid") (`var_sub') (`se_sub') (`ll_sub') (`ul_sub')
post `memhold' ("s2")        (`s2_step1') (`se_s2_step1') (`ll_s2_step1') (`ul_s2_step1')
postclose `memhold'

use "`results'", clear
export delimited using "`out_dir'/scen3_strata_coef_`sim_id'.csv", replace
