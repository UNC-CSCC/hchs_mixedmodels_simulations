version 19
set more off
set linesize 255

* Scenario 15: MI with EFFECTIVE sample size scaled weights (s1) as imputation-model covariates
* Usage: do scen15_mi.do <sim_id> <variant> <population>
*   variant 1: log(w_bg) + log(hh_s1)
*   variant 2: log(w_bg) + log(hh_s1) + log(sub_s1)
*   variant 3: log(w_bg) + log(hhsub_s1)   (hhsub = scaled HH x subject product weight)
*   population: dec2024 (v3data/Dec2024 samples; SEs clustered on block groups)
*               or oct2026 (v3data/Oct2026 samples; SEs clustered on households)
args sim_id variant population

if "`sim_id'" == "" {
    local sim_id "1"
}

if !inlist("`variant'", "1", "2", "3") {
    display as error "variant must be 1, 2, or 3"
    exit 198
}

if !inlist("`population'", "dec2024", "oct2026") {
    display as error "population must be dec2024 or oct2026"
    exit 198
}
local repo "/work/users/a/q/aquijano/mixed_sims/hchs_mixedmodels_simulations"
local data_dir "/work/users/b/e/beibo/v3data/Dec2024/sample_miss00pct"
local cluster_var "gl_bgid"
if "`population'" == "oct2026" {
    local data_dir "/work/users/b/e/beibo/v3data/Oct2026/sample_miss00pct"
    local cluster_var "gl_hhid"
}
local out_root "`repo'/data/estimates/`population'"

* Weight scaling method used for the imputation-model covariates
local s "s1"

* Ensure output directory exists
capture mkdir "`repo'/data/estimates"
capture mkdir "`out_root'"
capture mkdir "`out_root'/scen15_mi"
capture mkdir "`out_root'/scen15_mi/v`variant'"

import delimited "`data_dir'/samplemiss_`sim_id'.csv", clear varnames(1) bindquote(strict)
rename *, lower
destring _all, replace ignore("NA")

* Recode stratification variable to have values 1-4 and create indicator variables
gen int strat_recoded = mod(strat - 1, 4) + 1
gen byte strat_rec_2 = (strat_recoded == 2)
gen byte strat_rec_3 = (strat_recoded == 3)
gen byte strat_rec_4 = (strat_recoded == 4)
gen byte age_strat_new = (age_base >= 45)

* Recode IDs by group
egen long gl_bgid  = group(bgid)
egen long gl_hhid  = group(bgid hhid)
egen long gl_subid = group(bgid hhid subid)

* Filter valid design weights and subject identifiers
* (same analysis sample for all three variants; sub_s1 is always needed by the analysis model)
keep if !missing(w_sub, w_bghh, sub_s1, bgid, hhid, subid, x14, w_bg, hh_`s', sub_`s', hhsub_`s') ///
      & w_sub > 0 & w_bghh > 0 & sub_s1 > 0 & w_bg > 0 & hh_`s' > 0 & sub_`s' > 0 & hhsub_`s' > 0

* Log-transform weights for imputation
gen double log_w_bg  = log(w_bg)
gen double log_hh    = log(hh_`s')
gen double log_sub   = log(sub_`s')
gen double log_hhsub = log(hhsub_`s')

* Weight covariates entering the imputation model
if "`variant'" == "1" local wt_pred "log_w_bg log_hh"
if "`variant'" == "2" local wt_pred "log_w_bg log_hh log_sub"
if "`variant'" == "3" local wt_pred "log_w_bg log_hhsub"
display "Scenario 15 (`s'), variant `variant': imputation weight covariates = `wt_pred'"

* Mark unobserved visits under MAR missingness flag (miss_ind_mar == 1)
replace y_gfr = . if miss_ind_mar == 1
replace y_bmi = . if miss_ind_mar == 1

* Keep relevant variables
keep gl_bgid gl_hhid gl_subid bgid hhid subid strat_recoded strat_rec_2 strat_rec_3 strat_rec_4 age_strat_new x12 x14 x17 x18 w_sub w_bghh sub_s1 v_num y_gfr y_bmi x6 `wt_pred'

* ==============================================================================
* RESHAPE TO WIDE FORMAT FOR MULTIPLE IMPUTATION
* ==============================================================================
rename (y_gfr y_bmi x6) (y_gfr_ y_bmi_ x6_)
reshape wide y_gfr_ y_bmi_ x6_, i(gl_subid) j(v_num)

* Set up MI environment in wide format
mi set wide
mi register imputed y_bmi_2 y_gfr_2 y_bmi_3 y_gfr_3
mi register regular age_strat_new x12 x14 x17 x18 y_gfr_1 y_bmi_1 x6_2 x6_3 `wt_pred' strat_rec_2 strat_rec_3 strat_rec_4 w_sub w_bghh sub_s1 gl_hhid gl_bgid strat_recoded

* Imputation model: chained linear regression across wide visit variables (M = 5)
local M = 5
local mi_seed = 100000 + `sim_id'
mi impute chained ///
    (regress) y_bmi_2 y_gfr_2 y_bmi_3 y_gfr_3 = ///
    age_strat_new x12 x14 x17 x18 y_gfr_1 y_bmi_1 x6_2 x6_3 `wt_pred' strat_rec_2 strat_rec_3 strat_rec_4, ///
    add(`M') rseed(`mi_seed')

* Matrix to store parameter estimates (b_all) and sampling variances (v_all) across the M imputations
* Target terms: _cons, age_strat_new, x12, x17, x18, x6, y_bmi, var_hhid, var_subid, s2
local terms "_cons age_strat_new x12 x17 x18 x6 y_bmi var_hhid var_subid s2"
local n_terms : word count `terms'
matrix b_all = J(`n_terms', `M', 0)
matrix v_all = J(`n_terms', `M', 0)

* ==============================================================================
* RUN SCENARIO 9 TWO-STEP ESTIMATION ACROSS THE M IMPUTED DATASETS
* ==============================================================================
forvalues m = 1/`M' {
    preserve
        * Extract completed imputation m (in wide format)
        mi extract `m', clear

        * Reshape back to LONG format before running meglm
        reshape long y_gfr_ y_bmi_ x6_, i(gl_subid) j(v_num)
        rename (y_gfr_ y_bmi_ x6_) (y_gfr y_bmi x6)

        gen w_visit = 1

        * ----------------------------------------------------------------------
        * STEP 1: FIT WITH GLOBALLY SCALED WEIGHTS (mean = 1) -> fixed effects & s2
        * ----------------------------------------------------------------------
        egen double sum_w_sub_total = total(w_sub)
        gen  double w_sub_glob      = w_sub * (_N / sum_w_sub_total)

        capture meglm y_gfr age_strat_new x12 x17 x18 x6 y_bmi [pw=w_visit] || ///
            gl_hhid:  , pweight(w_bghh) || ///
            gl_subid: , pweight(w_sub_glob) ///
            vce(cluster `cluster_var') iterate(50) difficult
        * Stop rather than silently reuse estimates left over from a previous fit
        if _rc {
            display as error "Step 1 meglm failed (m = `m', rc = " _rc ")"
            exit _rc
        }
        display "FITCHECK sim=`sim_id' m=`m' step=1 converged=" e(converged)

        tempname rt1
        capture matrix `rt1' = r(table)

        foreach var in _cons age_strat_new x12 x17 x18 x6 y_bmi {
            local b_`var'  = _b[`var']
            local se_`var' = _se[`var']
            if missing(`se_`var'') local se_`var' = 0
            local v_`var'  = (`se_`var'')^2
        }

        local s2_step1    = .
        local se_s2_step1 = 0
        capture {
            local col_s2_step1 = colnumb(`rt1', "/:var(e.y_gfr)")
            local s2_step1     = `rt1'[1, `col_s2_step1']
            local se_s2_step1  = `rt1'[2, `col_s2_step1']
        }
        if missing(`se_s2_step1') local se_s2_step1 = 0
        local v_s2_step1   = (`se_s2_step1')^2

        * ----------------------------------------------------------------------
        * STEP 2: FIX FIXED EFFECTS AND RESIDUAL VARIANCE AT STEP 1 VALUES;
        *         ESTIMATE ONLY THE HOUSEHOLD AND SUBJECT VARIANCES
        *         Level 1 weight: sub_s1 (Pfeffermann effective size at subject level)
        *         Level 2 weight: w_bghh (unconditional household weight; keeps the block-group weight)
        * ----------------------------------------------------------------------
        constraint drop _all
        local c = 0
        foreach var in _cons age_strat_new x12 x17 x18 x6 y_bmi {
            local c = `c' + 1
            constraint `c' [y_gfr]`var' = `b_`var''
        }
        * Residual variance held at Step 1 so the reported s2 is the one used for var_hhid/var_subid
        local c = `c' + 1
        constraint `c' [/]var(e.y_gfr) = `s2_step1'

        capture meglm y_gfr age_strat_new x12 x17 x18 x6 y_bmi [pw=w_visit] || ///
            gl_hhid:  , pweight(w_bghh) || ///
            gl_subid: , pweight(sub_s1) ///
            vce(cluster `cluster_var') ///
            constraints(1/`c') iterate(50) difficult
        if _rc {
            display as error "Step 2 meglm failed (m = `m', rc = " _rc ")"
            exit _rc
        }
        display "FITCHECK sim=`sim_id' m=`m' step=2 converged=" e(converged)

        tempname rt2
        capture matrix `rt2' = r(table)

        local var_hh = .
        local se_hh  = 0
        capture {
            local col_hh  = colnumb(`rt2', "/:var(_cons[gl_hhid])")
            local var_hh  = `rt2'[1, `col_hh']
            local se_hh   = `rt2'[2, `col_hh']
        }
        if missing(`se_hh') local se_hh = 0
        local v_hh    = (`se_hh')^2

        local var_sub = .
        local se_sub  = 0
        capture {
            local col_sub = colnumb(`rt2', "/:var(_cons[gl_hhid>gl_subid])")
            local var_sub = `rt2'[1, `col_sub']
            local se_sub  = `rt2'[2, `col_sub']
        }
        if missing(`se_sub') local se_sub = 0
        local v_sub   = (`se_sub')^2

        * Store into matrices
        matrix b_all[1, `m'] = `b__cons'
        matrix v_all[1, `m'] = `v__cons'

        matrix b_all[2, `m'] = `b_age_strat_new'
        matrix v_all[2, `m'] = `v_age_strat_new'

        matrix b_all[3, `m'] = `b_x12'
        matrix v_all[3, `m'] = `v_x12'

        matrix b_all[4, `m'] = `b_x17'
        matrix v_all[4, `m'] = `v_x17'

        matrix b_all[5, `m'] = `b_x18'
        matrix v_all[5, `m'] = `v_x18'

        matrix b_all[6, `m'] = `b_x6'
        matrix v_all[6, `m'] = `v_x6'

        matrix b_all[7, `m'] = `b_y_bmi'
        matrix v_all[7, `m'] = `v_y_bmi'

        matrix b_all[8, `m'] = `var_hh'
        matrix v_all[8, `m'] = `v_hh'

        matrix b_all[9, `m'] = `var_sub'
        matrix v_all[9, `m'] = `v_sub'

        matrix b_all[10, `m'] = `s2_step1'
        matrix v_all[10, `m'] = `v_s2_step1'
    restore
}

* ==============================================================================
* POOL RESULTS USING RUBIN'S RULES ACROSS THE M IMPUTATIONS
* ==============================================================================
tempname memhold
tempfile results
postfile `memhold' str32 term double estimate double std_error double conf_low double conf_high ///
    using "`results'", replace

local z975 = invnormal(0.975)

forvalues i = 1/`n_terms' {
    local t_name : word `i' of `terms'
    local out_name = "`t_name'"
    if "`t_name'" == "_cons" local out_name "(Intercept)"

    * Mean point estimate (Q_bar)
    local q_sum = 0
    local u_sum = 0
    forvalues m = 1/`M' {
        local q_sum = `q_sum' + b_all[`i', `m']
        local u_sum = `u_sum' + v_all[`i', `m']
    }
    local q_bar = `q_sum' / `M'
    local u_bar = `u_sum' / `M'

    * Between-imputation variance (B)
    local b_sum = 0
    forvalues m = 1/`M' {
        local diff = b_all[`i', `m'] - `q_bar'
        local b_sum = `b_sum' + (`diff')^2
    }
    local B = `b_sum' / (`M' - 1)

    * Total variance (T) and standard error
    local T = `u_bar' + (1 + 1 / `M') * `B'
    local se_pooled = sqrt(`T')
    local ci_low  = `q_bar' - `z975' * `se_pooled'
    local ci_high = `q_bar' + `z975' * `se_pooled'

    * Variance components (terms 8-10): interval pooled on the log scale and back-transformed,
    * matching how meglm builds variance intervals in the complete-data scen9 (log SE = SE / estimate)
    if `i' >= 8 {
        local ql_sum = 0
        local ul_sum = 0
        forvalues m = 1/`M' {
            local ql_sum = `ql_sum' + ln(b_all[`i', `m'])
            local ul_sum = `ul_sum' + v_all[`i', `m'] / (b_all[`i', `m'])^2
        }
        local ql_bar = `ql_sum' / `M'
        local ul_bar = `ul_sum' / `M'
        local bl_sum = 0
        forvalues m = 1/`M' {
            local bl_sum = `bl_sum' + (ln(b_all[`i', `m']) - `ql_bar')^2
        }
        local T_log = `ul_bar' + (1 + 1 / `M') * (`bl_sum' / (`M' - 1))
        local ci_low  = exp(`ql_bar' - `z975' * sqrt(`T_log'))
        local ci_high = exp(`ql_bar' + `z975' * sqrt(`T_log'))
    }

    post `memhold' ("`out_name'") (`q_bar') (`se_pooled') (`ci_low') (`ci_high')
}

postclose `memhold'

use "`results'", clear
export delimited using "`out_root'/scen15_mi/v`variant'/scen15_mi_coef_`sim_id'.csv", replace

* Per-imputation estimates and sampling variances, so the pooling can be redone without refitting
capture mkdir "`out_root'/scen15_mi/v`variant'/per_imputation"
tempname implog
tempfile imp_results
postfile `implog' str32 term int imputation double estimate double variance using "`imp_results'", replace
forvalues i = 1/`n_terms' {
    local t_name : word `i' of `terms'
    if "`t_name'" == "_cons" local t_name "(Intercept)"
    forvalues m = 1/`M' {
        post `implog' ("`t_name'") (`m') (b_all[`i', `m']) (v_all[`i', `m'])
    }
}
postclose `implog'
use "`imp_results'", clear
export delimited using "`out_root'/scen15_mi/v`variant'/per_imputation/scen15_mi_imp_`sim_id'.csv", replace
