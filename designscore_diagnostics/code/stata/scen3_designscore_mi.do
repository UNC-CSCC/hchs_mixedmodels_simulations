version 19
set more off
set linesize 255
* Paired MI comparison: scenario 4 s1/v2 imputations, original scen3 vs design score.
* Called by run_designscore_mi.py. All input paths and outputs are explicit.
args sim_id diagnostic_root run_label M
cd `"`diagnostic_root'"'
do code/stata/mi_rename_compat.do

capture program drop _dsmi_run
program define _dsmi_run
    version 19
    args sim_id run_label M
    local sim = real("`sim_id'")
    local M = real("`M'")
    if missing(`sim',`M') | `sim'<1 | `M'<2 exit 198
    local out "data/estimates/oct2026/`run_label'"
    local sample "../../data/Oct2026/sample_miss00pct/samplemiss_`sim'.csv"
    local impfile "`out'/imputations/mi_`sim'.dta"
    local coeffile "`out'/pooled/mi_coef_`sim'.csv"
    capture confirm file "`coeffile'"
    if !_rc exit 602
    log using "logs/oct2026/`run_label'/mi_fit_`sim'.log", text name(dsmi) replace

    * Exact subject selection probabilities; unchanged across imputations.
    import delimited "data/design_frames/oct2026_subject_counts.csv", clear varnames(1) case(lower)
    keep if population == "oct2026" & sim_id == `sim'
    assert _N == 1 & ids_verified == "TRUE"
    assert eligible_younger >= selected_younger & selected_younger >= 2
    assert eligible_older >= selected_older & selected_older >= 2
    tempname py po pyy poo probs
    scalar `py' = selected_younger[1]/eligible_younger[1]
    scalar `po' = selected_older[1]/eligible_older[1]
    scalar `pyy' = selected_younger[1]*(selected_younger[1]-1)/(eligible_younger[1]*(eligible_younger[1]-1))
    scalar `poo' = selected_older[1]*(selected_older[1]-1)/(eligible_older[1]*(eligible_older[1]-1))
    matrix `probs' = (`py',`po',`pyy',`poo',`py'*`po')
    local selected_y = selected_younger[1]
    local selected_o = selected_older[1]

    * Reuse the exact saved imputation set if a fit must be resumed.
    capture confirm file "`impfile'"
    if _rc {
        import delimited "data/design_frames/oct2026_households.csv", clear varnames(1) case(lower)
        assert population == "oct2026"
        keep bgid hhid strat m_h p_h
        isid bgid hhid
        rename strat frame_strat
        tempfile frame
        save "`frame'"
        import delimited "`sample'", clear varnames(1) case(lower) bindquote(strict)
        assert !missing(y_gfr,y_bmi,x6,x12,x14,x17,x18,age_base,w_sub,w_bghh,sub_s1,hh_s1,w_bg)
        assert inlist(miss_ind_mar,0,1)
        assert miss_ind_mar == 0 if v_num == 1
        assert base_miss == 0
        assert w_sub>0 & w_bghh>0 & sub_s1>0 & hh_s1>0 & w_bg>0
        isid bgid hhid subid v_num
        bysort bgid hhid subid (v_num): assert _N==3 & v_num==_n
        gen byte age_strat_new = age_base >= 45
        quietly count if v_num==1 & age_strat_new==0
        assert r(N)==`selected_y'
        quietly count if v_num==1 & age_strat_new==1
        assert r(N)==`selected_o'
        merge m:1 bgid hhid using "`frame'", keep(master match)
        assert _merge==3 & strat==frame_strat
        drop _merge frame_strat
        assert m_h>=1 & p_h>0 & p_h<=1
        egen long gl_hhid = group(bgid hhid)
        egen long gl_subid = group(bgid hhid subid)
        egen long gl_psu = group(strat bgid)
        gen int strat_recoded = mod(strat-1,4)+1
        forvalues k=2/4 {
            gen byte strat_rec_`k' = strat_recoded==`k'
        }
        gen double log_w_bg = log(w_bg)
        gen double log_hh = log(hh_s1)
        gen double log_sub = log(sub_s1)
        * No numeric outcome from a flagged visit enters the imputation fit.
        replace y_gfr = . if miss_ind_mar==1
        replace y_bmi = . if miss_ind_mar==1
        keep bgid hhid subid strat gl_hhid gl_subid gl_psu m_h p_h age_strat_new x12 x14 x17 x18 ///
            w_sub w_bghh sub_s1 v_num y_gfr y_bmi x6 log_w_bg log_hh log_sub strat_rec_2 strat_rec_3 strat_rec_4 miss_ind_mar
        rename y_gfr y_gfr_
        rename y_bmi y_bmi_
        rename x6 x6_
        rename miss_ind_mar mar_
        reshape wide y_gfr_ y_bmi_ x6_ mar_, i(gl_subid) j(v_num)
        sort gl_subid
        mi set wide
        mi register imputed y_bmi_2 y_gfr_2 y_bmi_3 y_gfr_3
        mi register regular bgid hhid subid strat gl_hhid gl_subid gl_psu m_h p_h age_strat_new x12 x14 x17 x18 ///
            w_sub w_bghh sub_s1 y_gfr_1 y_bmi_1 x6_1 x6_2 x6_3 log_w_bg log_hh log_sub strat_rec_2 strat_rec_3 strat_rec_4 mar_1 mar_2 mar_3
        local seed = 100000+`sim'
        mi impute chained (regress) y_bmi_2 y_gfr_2 y_bmi_3 y_gfr_3 = ///
            age_strat_new x12 x14 x17 x18 y_gfr_1 y_bmi_1 x6_2 x6_3 ///
            log_w_bg log_hh log_sub strat_rec_2 strat_rec_3 strat_rec_4, add(`M') rseed(`seed')
        save "`impfile'"
    }
    else use "`impfile'", clear
    quietly mi describe
    assert r(M)==`M'
    tempfile imputed
    save "`imputed'"
    quietly do code/stata/variance_designscore.mata
    local terms "(Intercept) age_strat_new x12 x17 x18 x6 y_bmi var_hhid var_subid s2"
    tempname per auditpost
    tempfile perdata auditdata
    postfile `per' long sim_id int imputation str12 method str32 term double estimate std_error ///
        str100 input_signature byte inference_usable using "`perdata'", replace
    postfile `auditpost' long sim_id int imputation str100 input_signature long n_subjects n_households n_psu ///
        double missing_v2 missing_v3 var_hhid_step1 var_subid_step1 s2_step1 ///
        boundary_mode se_status kkt_error step1_score_error jacobian_condition original_converged ///
        using "`auditdata'", replace

    forvalues m=1/`M' {
        use "`imputed'", clear
        mi extract `m', clear
        reshape long y_gfr_ y_bmi_ x6_ mar_, i(gl_subid) j(v_num)
        rename y_gfr_ y_gfr
        rename y_bmi_ y_bmi
        rename x6_ x6
        rename mar_ miss_ind_mar
        assert !missing(y_gfr,y_bmi,x6,x12,x14,x17,x18)
        bysort gl_subid (v_num): assert _N==3 & v_num==_n
        gen double ds_one = 1
        gen double w_visit = 1
        egen double sum_w_sub_total = total(w_sub)
        tempname scale rt1 initial dsest dsaudit orig
        scalar `scale' = _N/sum_w_sub_total[1]
        gen double w_sub_glob = w_sub*`scale'
        sort gl_hhid gl_subid v_num
        quietly _datasignature y_gfr ds_one age_strat_new x12 x17 x18 x6 y_bmi gl_hhid gl_subid gl_psu strat m_h p_h w_sub w_bghh sub_s1
        local signature "`r(datasignature)'"
        local nsub = _N/3
        quietly summarize miss_ind_mar if v_num==2, meanonly
        local miss2 = r(mean)
        quietly summarize miss_ind_mar if v_num==3, meanonly
        local miss3 = r(mean)

        * A shared first-stage fit is valid because both methods use this same fit.
        capture noisily meglm y_gfr age_strat_new x12 x17 x18 x6 y_bmi [pw=w_visit] || ///
            gl_hhid: , pweight(w_bghh) || gl_subid: , pweight(w_sub_glob) vce(cluster gl_hhid) iterate(300) nolog
        if _rc exit _rc
        if e(converged)!=1 exit 430
        matrix `rt1' = r(table)
        assert e(N)==_N
        matrix `initial' = J(1,11,.)
        matrix `orig' = J(10,2,.)
        local j=0
        foreach v in _cons age_strat_new x12 x17 x18 x6 y_bmi {
            local ++j
            matrix `initial'[1,`j'] = _b[`v']
            matrix `orig'[`j',1] = _b[`v']
            matrix `orig'[`j',2] = _se[`v']
        }
        local hhcol = colnumb(`rt1',"/:var(_cons[gl_hhid])")
        local subcol = colnumb(`rt1',"/:var(_cons[gl_hhid>gl_subid])")
        local ecol = colnumb(`rt1',"/:var(e.y_gfr)")
        if missing(`hhcol',`subcol',`ecol') exit 498
        matrix `initial'[1,8] = `rt1'[1,`hhcol']
        matrix `initial'[1,9] = `rt1'[1,`subcol']
        matrix `initial'[1,10] = `rt1'[1,`ecol']
        matrix `initial'[1,11] = `scale'
        matrix `orig'[10,1] = `rt1'[1,`ecol']
        matrix `orig'[10,2] = `rt1'[2,`ecol']

        mata: ds_run("y_gfr ds_one age_strat_new x12 x17 x18 x6 y_bmi", ///
            "gl_hhid gl_psu strat m_h w_bghh p_h w_sub age_strat_new", ///
            "`initial'","`probs'","`dsest'","`dsaudit'")
        quietly _datasignature y_gfr ds_one age_strat_new x12 x17 x18 x6 y_bmi gl_hhid gl_subid gl_psu strat m_h p_h w_sub w_bghh sub_s1
        assert "`r(datasignature)'"=="`signature'"

        * Original scenario 3: constrained second-stage meglm, same completed data.
        constraint drop _all
        local j=0
        foreach v in _cons age_strat_new x12 x17 x18 x6 y_bmi {
            local ++j
            local value : display %21.17g `initial'[1,`j']
            constraint `j' [y_gfr]`v' = `value'
        }
        local value : display %21.17g `initial'[1,10]
        constraint 8 [/]var(e.y_gfr) = `value'
        capture noisily meglm y_gfr age_strat_new x12 x17 x18 x6 y_bmi [pw=w_visit] || ///
            gl_hhid: , pweight(w_bghh) || gl_subid: , pweight(sub_s1) ///
            vce(cluster gl_hhid) constraints(1/8) iterate(300) nolog
        if _rc exit _rc
        if e(converged)!=1 exit 430
        tempname rt2
        matrix `rt2' = r(table)
        local hhcol = colnumb(`rt2',"/:var(_cons[gl_hhid])")
        local subcol = colnumb(`rt2',"/:var(_cons[gl_hhid>gl_subid])")
        matrix `orig'[8,1] = `rt2'[1,`hhcol']
        matrix `orig'[8,2] = `rt2'[2,`hhcol']
        matrix `orig'[9,1] = `rt2'[1,`subcol']
        matrix `orig'[9,2] = `rt2'[2,`subcol']
        quietly _datasignature y_gfr ds_one age_strat_new x12 x17 x18 x6 y_bmi gl_hhid gl_subid gl_psu strat m_h p_h w_sub w_bghh sub_s1
        assert "`r(datasignature)'"=="`signature'"
        forvalues j=1/10 {
            local term : word `j' of `terms'
            local valid = !missing(`orig'[`j',1],`orig'[`j',2]) & `orig'[`j',2]>0
            if `j'>=8 local valid = `valid' & `orig'[`j',1]>0
            post `per' (`sim') (`m') ("original") ("`term'") (`orig'[`j',1]) (`orig'[`j',2]) ("`signature'") (`valid')
            local valid = `dsaudit'[1,6]==0 & !missing(`dsest'[`j',1],`dsest'[`j',2]) & `dsest'[`j',2]>0
            if `j'>=8 local valid = `valid' & `dsest'[`j',1]>0
            post `per' (`sim') (`m') ("designscore") ("`term'") (`dsest'[`j',1]) (`dsest'[`j',2]) ("`signature'") (`valid')
        }
        post `auditpost' (`sim') (`m') ("`signature'") (`nsub') (`dsaudit'[1,10]) (`dsaudit'[1,9]) ///
            (`miss2') (`miss3') (`initial'[1,8]) (`initial'[1,9]) (`initial'[1,10]) ///
            (`dsaudit'[1,1]) (`dsaudit'[1,6]) (`dsaudit'[1,2]) (`dsaudit'[1,7]) (`dsaudit'[1,8]) (1)
        display "PAIRED_MI_FIT sim=`sim' imputation=`m' signature=`signature'"
    }
    postclose `per'
    postclose `auditpost'
    use "`perdata'", clear
    format estimate std_error %24.16g
    export delimited using "`out'/per_imputation/mi_estimates_`sim'.csv", replace
    use "`auditdata'", clear
    format missing_v2-jacobian_condition %24.16g
    export delimited using "`out'/audits/mi_audit_`sim'.csv", replace
    * Pooling is performed independently in Python, with all M required.
    display "PAIRED_MI_COMPLETE sim=`sim' M=`M'"
end

capture noisily _dsmi_run "`sim_id'" "`run_label'" "`M'"
local rc=_rc
capture log close dsmi
capture program drop _dsmi_run
if `rc' exit `rc'
