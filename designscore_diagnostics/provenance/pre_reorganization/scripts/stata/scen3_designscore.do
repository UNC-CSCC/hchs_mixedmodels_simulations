version 19
set more off
set linesize 255

* New variant. The original scen3.do is not called or modified.
* Usage (run from the submodule root, or supply repo_root):
* do scripts/stata/scen3_designscore.do <sim_id> <population> <exact|bernoulli> ///
*     [repo_root] [sample_directory] [household_frame_csv] [subject_counts_csv]
* See SCEN3_DESIGNSCORE.md for the required metadata and inference limitations.
args sim_id population probability_mode repo_root sample_directory frame_file counts_file

capture program drop _scen3_designscore
program define _scen3_designscore
    version 19
    args sim_id population probability_mode repo_root sample_directory frame_file counts_file

    if !inlist("`population'", "dec2024", "oct2026") {
        display as error "population must be dec2024 or oct2026"
        exit 198
    }
    if !inlist("`probability_mode'", "exact", "bernoulli") {
        display as error "Explicitly choose exact or bernoulli probability mode."
        exit 198
    }
    local sim = real("`sim_id'")
    if missing(`sim') | `sim' < 1 | `sim' != floor(`sim') {
        display as error "sim_id must be a positive integer"
        exit 198
    }
    local sim_id = strtrim(string(`sim', "%12.0f"))
    if `"`repo_root'"' == "" local repo_root `"`c(pwd)'"'
    local release = cond("`population'" == "oct2026", "Oct2026", "Dec2024")
    if `"`sample_directory'"' == "" {
        local sample_directory `"`repo_root'/../data/`release'/sample_miss00pct"'
    }
    if `"`frame_file'"' == "" {
        local frame_file `"`repo_root'/data/design_frames/`population'_households.csv"'
    }
    if `"`counts_file'"' == "" {
        local counts_file `"`repo_root'/data/design_frames/`population'_subject_counts.csv"'
    }
    local sample_file `"`sample_directory'/samplemiss_`sim_id'.csv"'
    local solver_file `"`repo_root'/scripts/stata/variance_designscore.mata"'
    confirm file `"`sample_file'"'
    confirm file `"`frame_file'"'
    confirm file `"`solver_file'"'
    if "`probability_mode'" == "exact" confirm file `"`counts_file'"'

    local out_dir `"`repo_root'/data/estimates/`population'/scen3_designscore/`probability_mode'"'
    capture mkdir `"`repo_root'/data/estimates"'
    capture mkdir `"`repo_root'/data/estimates/`population'"'
    capture mkdir `"`repo_root'/data/estimates/`population'/scen3_designscore"'
    capture mkdir `"`out_dir'"'
    local coefficient_file `"`out_dir'/scen3_designscore_coef_`sim_id'.csv"'
    local diagnostic_file `"`out_dir'/scen3_designscore_audit_`sim_id'.csv"'
    foreach path in coefficient_file diagnostic_file {
        capture confirm file `"``path''"'
        if !_rc {
            display as error "Output already exists: ``path''. Choose a fresh output location."
            exit 602
        }
    }
    capture log close dsdesignscore
    log using `"`out_dir'/scen3_designscore_`sim_id'.log"', text name(dsdesignscore) replace
    display "DESIGNSCORE sample: `sample_file'"
    display "DESIGNSCORE frame: `frame_file'"
    display "DESIGNSCORE probability_mode=`probability_mode'"

    * Import and validate the full-household lookup BEFORE fitting any model.
    import delimited `"`frame_file'"', clear varnames(1) case(lower) bindquote(strict)
    if _N == 0 {
        display as error "The household frame is empty."
        exit 459
    }
    confirm string variable population
    assert population == "`population'"
    confirm numeric variable bgid hhid strat m_h p_h
    assert !missing(bgid, hhid, strat, m_h, p_h)
    assert bgid >= 1 & bgid == floor(bgid) & hhid >= 1 & hhid == floor(hhid)
    assert m_h >= 1 & m_h == floor(m_h)
    assert p_h > 0 & p_h <= 1
    assert inrange(strat, 1, 8) & strat == floor(strat)
    isid bgid hhid
    keep bgid hhid strat m_h p_h
    gen byte frame_strat = strat
    drop strat
    tempfile household_frame
    save "`household_frame'"

    tempname py po pyy poo pyo eligible_y eligible_o selected_y selected_o
    scalar `eligible_y' = .
    scalar `eligible_o' = .
    scalar `selected_y' = .
    scalar `selected_o' = .
    if "`probability_mode'" == "exact" {
        import delimited `"`counts_file'"', clear varnames(1) case(lower) bindquote(strict)
        confirm string variable population
        confirm numeric variable sim_id eligible_younger eligible_older selected_younger selected_older
        isid population sim_id
        keep if population == "`population'" & sim_id == `sim'
        if _N != 1 {
            display as error "Exactly one subject-count row is required for this population and replicate."
            exit 459
        }
        foreach g in younger older {
            assert !missing(eligible_`g', selected_`g')
            assert eligible_`g' >= 2 & eligible_`g' == floor(eligible_`g')
            assert selected_`g' >= 2 & selected_`g' == floor(selected_`g')
            assert selected_`g' <= eligible_`g'
        }
        scalar `eligible_y' = eligible_younger[1]
        scalar `eligible_o' = eligible_older[1]
        scalar `selected_y' = selected_younger[1]
        scalar `selected_o' = selected_older[1]
        scalar `py' = `selected_y' / `eligible_y'
        scalar `po' = `selected_o' / `eligible_o'
        scalar `pyy' = `selected_y' * (`selected_y' - 1) / (`eligible_y' * (`eligible_y' - 1))
        scalar `poo' = `selected_o' * (`selected_o' - 1) / (`eligible_o' * (`eligible_o' - 1))
    }
    else {
        scalar `py' = 0.3575
        scalar `po' = 0.55
        scalar `pyy' = `py'^2
        scalar `poo' = `po'^2
        display "DESIGNSCORE: nominal Bernoulli subject probabilities; pair probabilities are approximate."
    }
    scalar `pyo' = `py' * `po'

    import delimited `"`sample_file'"', clear varnames(1) case(lower) bindquote(strict)
    if _N == 0 {
        display as error "The sample is empty."
        exit 459
    }
    confirm numeric variable strat bgid hhid subid v_num age_base x12 x17 x18 x6 x14 y_bmi y_gfr w_sub w_bghh sub_s1
    assert !missing(strat, bgid, hhid, subid, v_num, age_base, x12, x17, x18, x6, x14, y_bmi, y_gfr, w_sub, w_bghh, sub_s1)
    assert w_sub > 0 & w_bghh > 0 & sub_s1 > 0
    assert inrange(strat,1,8) & strat == floor(strat)
    assert bgid >= 1 & hhid >= 1 & subid >= 1
    assert bgid == floor(bgid) & hhid == floor(hhid) & subid == floor(subid)
    assert inrange(v_num,1,3) & v_num == floor(v_num)
    capture confirm variable base_miss
    if !_rc assert base_miss == 0
    capture confirm variable dat_num
    if !_rc assert dat_num == `sim'
    isid bgid hhid subid v_num
    bysort bgid hhid subid (v_num): assert _N == 3 & v_num == _n
    foreach v in age_base x12 x17 x18 x14 w_sub w_bghh sub_s1 strat {
        by bgid hhid subid: assert `v' == `v'[1]
    }
    bysort bgid hhid: assert w_bghh == w_bghh[1] & strat == strat[1]
    bysort bgid: assert strat == strat[1]
    gen byte age_strat_new = age_base >= 45
    assert abs(w_sub - cond(age_strat_new, 1/0.55, 1/0.3575)) < 1e-7
    if "`probability_mode'" == "exact" {
        quietly count if v_num == 1 & age_strat_new == 0
        assert r(N) == `selected_y'
        quietly count if v_num == 1 & age_strat_new == 1
        assert r(N) == `selected_o'
    }

    merge m:1 bgid hhid using "`household_frame'", keep(master match)
    assert _merge == 3
    assert strat == frame_strat
    drop _merge frame_strat
    bysort bgid hhid: assert _N/3 <= m_h
    egen long gl_hhid = group(bgid hhid)
    egen long gl_subid = group(bgid hhid subid)
    egen long gl_psu = group(strat bgid)
    gen double ds_one = 1
    gen double w_visit = 1
    egen double sum_w_sub_total = total(w_sub)
    gen double w_sub_glob = w_sub * (_N / sum_w_sub_total)
    tempname global_scale
    scalar `global_scale' = _N / sum_w_sub_total[1]
    sort gl_hhid gl_subid v_num

    * Define the helper if this is a fresh Stata session. Nothing fits on load.
    mata: st_local("ds_loaded", strofreal(findexternal("ds_run()") != NULL))
    if !`ds_loaded' quietly do `"`solver_file'"'

    * STEP 1: same model and weights as scen3.do; complete-data benchmark.
    * Numeric values at flagged missed visits are intentionally retained here.
    meglm y_gfr age_strat_new x12 x17 x18 x6 y_bmi [pw=w_visit] || ///
        gl_hhid: , pweight(w_bghh) || ///
        gl_subid: , pweight(w_sub_glob) vce(cluster gl_hhid)
    tempname rt1 initial probabilities estimates audit
    matrix `rt1' = r(table)
    if e(converged) != 1 {
        display as error "Step 1 did not converge; no coefficient file will be written."
        exit 430
    }
    assert e(N) == _N
    matrix `initial' = J(1,11,.)
    local j = 0
    foreach v in _cons age_strat_new x12 x17 x18 x6 y_bmi {
        local ++j
        matrix `initial'[1,`j'] = _b[`v']
    }
    local hhcol = colnumb(`rt1', "/:var(_cons[gl_hhid])")
    local subcol = colnumb(`rt1', "/:var(_cons[gl_hhid>gl_subid])")
    local ecol = colnumb(`rt1', "/:var(e.y_gfr)")
    if missing(`hhcol', `subcol', `ecol') {
        display as error "Expected Step 1 variance columns are absent from r(table)."
        exit 498
    }
    matrix `initial'[1,8] = `rt1'[1,`hhcol']
    matrix `initial'[1,9] = `rt1'[1,`subcol']
    matrix `initial'[1,10] = `rt1'[1,`ecol']
    matrix `initial'[1,11] = `global_scale'
    matrix `probabilities' = (`py', `po', `pyy', `poo', `pyo')

    * STEP 2: weighted census scores, positive interior and explicit boundaries.
    mata: ds_run("y_gfr ds_one age_strat_new x12 x17 x18 x6 y_bmi", ///
        "gl_hhid gl_psu strat m_h w_bghh p_h w_sub age_strat_new", ///
        "`initial'", "`probabilities'", "`estimates'", "`audit'")
    matrix colnames `estimates' = estimate std_error conf_low conf_high
    matrix rownames `estimates' = intercept age_strat_new x12 x17 x18 x6 y_bmi var_hhid var_subid s2
    matrix colnames `audit' = boundary_mode kkt_error objective attempts accepted se_status step1_score_error jacobian_condition n_psu n_households n_subjects score_hh score_sub
    matrix list `audit'

    local se_status = `audit'[1,6]
    local boundary_mode = `audit'[1,1]
    local fit_status "OK"
    if `se_status' == 6 local fit_status "BOUNDARY_POINT_ONLY"
    else if `se_status' != 0 local fit_status "SE_FAILED_POINT_ONLY"
    if `se_status' != 0 {
        display as error "`fit_status': estimates retained; all SEs/CIs are missing. See audit and documentation."
    }

    * Keep legacy coefficient column names, adding status and method provenance.
    tempname post
    tempfile coefficient_rows coefficient_csv diagnostic_csv
    postfile `post' str32 term double estimate double std_error double conf_low double conf_high using "`coefficient_rows'"
    local terms "(Intercept) age_strat_new x12 x17 x18 x6 y_bmi var_hhid var_subid s2"
    forvalues j = 1/10 {
        local term : word `j' of `terms'
        post `post' ("`term'") (`estimates'[`j',1]) (`estimates'[`j',2]) (`estimates'[`j',3]) (`estimates'[`j',4])
    }
    postclose `post'
    use "`coefficient_rows'", clear
    gen long sim_id = `sim'
    gen str8 population = "`population'"
    gen str9 probability_mode = "`probability_mode'"
    gen str32 fit_status = "`fit_status'"
    gen str32 se_method = "stacked_stratified_psu_wr"
    gen byte se_status = `se_status'
    gen byte boundary_mode = `boundary_mode'
    gen byte inference_usable = se_status == 0
    format estimate std_error conf_low conf_high %24.16g
    export delimited using "`coefficient_csv'", replace

    clear
    quietly set obs 1
    svmat double `audit', names(col)
    gen long sim_id = `sim'
    gen str8 population = "`population'"
    gen str9 probability_mode = "`probability_mode'"
    gen str32 fit_status = "`fit_status'"
    gen str32 se_method = "stacked_stratified_psu_wr"
    gen double p_younger = `py'
    gen double p_older = `po'
    gen double pair_younger = `pyy'
    gen double pair_older = `poo'
    gen double pair_mixed = `pyo'
    gen double eligible_younger = `eligible_y'
    gen double eligible_older = `eligible_o'
    gen double step1_var_hhid = `initial'[1,8]
    gen double step1_var_subid = `initial'[1,9]
    gen double step1_s2 = `initial'[1,10]
    gen double global_weight_scale = `global_scale'
    quietly ds, has(type numeric)
    format `r(varlist)' %24.16g
    export delimited using "`diagnostic_csv'", replace
    * copy without replace: previously saved outputs are never overwritten.
    copy "`diagnostic_csv'" `"`diagnostic_file'"'
    copy "`coefficient_csv'" `"`coefficient_file'"'
    display "DESIGNSCORE `fit_status': `coefficient_file'"
end

capture noisily _scen3_designscore "`sim_id'" "`population'" "`probability_mode'" ///
    `"`repo_root'"' `"`sample_directory'"' `"`frame_file'"' `"`counts_file'"'
local rc = _rc
capture log close dsdesignscore
capture program drop _scen3_designscore
if `rc' {
    display as error "scen3_designscore failed with return code `rc'. See the log; no successful fit is implied."
    exit `rc'
}
