version 19
set more off
set linesize 255
* Run ONLY through run_designscore_mlmi.py. This file generates imputations;
* it never calls either analysis estimator. The Python driver checks chains
* before promoting the candidate MI set and calling the unchanged paired fit.
args sim_id diagnostic_root label M chains burn between every ridge hh_prior psu_prior person_prior
cd `"`diagnostic_root'"'
do code/stata/mi_rename_compat.do

capture program drop _mlmi_generate
program define _mlmi_generate
    version 19
    args sim_id label M chains burn between every ridge hh_prior psu_prior person_prior
    local sim=real("`sim_id'")
    if !regexm("`label'","^mlmi_[A-Za-z0-9_]+$") exit 198
    foreach n in sim M chains burn between every {
        if missing(``n'') | ``n''!=floor(``n'') | ``n''<1 exit 198
    }
    if `chains'<2 | mod(`M',`chains')!=0 exit 198
    if `M'/`chains'*`between'/`every'<100 exit 198
    if mod(`burn',`every') | mod(`between',`every') exit 198
    if min(`ridge',`hh_prior',`psu_prior',`person_prior')<=0 exit 198
    local out "data/estimates/oct2026/`label'"
    local candidate "`out'/imputations/candidate_`sim'.dta"
    confirm new file "`candidate'"
    confirm new file "`out'/imputations/mi_`sim'.dta"
    log using "logs/oct2026/`label'/mlmi_generate_`sim'.log", text name(mlmi)

    import delimited "data/design_frames/oct2026_households.csv", clear varnames(1) case(lower)
    assert population=="oct2026"
    keep bgid hhid strat m_h p_h
    isid bgid hhid
    rename strat frame_strat
    tempfile frame base completed
    save "`frame'"
    import delimited "../../data/Oct2026/sample_miss00pct/samplemiss_`sim'.csv", ///
        clear varnames(1) case(lower) bindquote(strict) asdouble
    assert !missing(y_gfr,y_bmi,x6,x12,x14,x17,x18,age_base,w_sub,w_bghh,sub_s1,hh_s1,w_bg)
    assert base_miss==0 & inlist(miss_ind_mar,0,1)
    assert miss_ind_mar==0 if v_num==1
    assert w_sub>0 & w_bghh>0 & sub_s1>0 & hh_s1>0 & w_bg>0
    isid bgid hhid subid v_num
    bysort bgid hhid subid (v_num): assert _N==3 & v_num==_n
    gen byte age_strat_new=age_base>=45
    merge m:1 bgid hhid using "`frame'", keep(master match)
    assert _merge==3 & strat==frame_strat
    drop _merge frame_strat
    egen long gl_hhid=group(bgid hhid)
    egen long gl_subid=group(bgid hhid subid)
    egen long gl_psu=group(strat bgid)
    gen int strat_recoded=mod(strat-1,4)+1
    forvalues k=2/4 {
        gen byte strat_rec_`k'=strat_recoded==`k'
    }
    forvalues k=2/8 {
        gen byte mlmi_strat_`k'=strat==`k'
    }
    gen double log_w_bg=log(w_bg)
    gen double log_hh=log(hh_s1)
    gen double log_sub=log(sub_s1)
    gen double mlmi_x17_age=x17*age_strat_new
    forvalues k=1/3 {
        gen double mlmi_x12_strat`k'=x12*(strat_recoded==`k')
    }
    * Mask before scaling, initializing or calculating any sampler statistic.
    replace y_gfr=. if miss_ind_mar==1
    replace y_bmi=. if miss_ind_mar==1
    keep bgid hhid subid strat gl_hhid gl_subid gl_psu m_h p_h age_base age_strat_new ///
        x12 x14 x17 x18 w_sub w_bghh sub_s1 v_num y_gfr y_bmi x6 ///
        log_w_bg log_hh log_sub strat_rec_2 strat_rec_3 strat_rec_4 miss_ind_mar mlmi_*
    rename y_gfr y_gfr_
    rename y_bmi y_bmi_
    rename x6 x6_
    rename miss_ind_mar mar_
    reshape wide y_gfr_ y_bmi_ x6_ mar_, i(gl_subid) j(v_num)
    sort gl_psu gl_hhid gl_subid
    recast double y_gfr_* y_bmi_*
    assert missing(y_gfr_2)==mar_2 & missing(y_bmi_2)==mar_2
    assert missing(y_gfr_3)==mar_3 & missing(y_bmi_3)==mar_3
    local ys "y_gfr_1 y_bmi_1 y_gfr_2 y_bmi_2 y_gfr_3 y_bmi_3"
    local xs "age_base age_strat_new x12 x14 x17 x18 x6_2 x6_3 mlmi_strat_2 mlmi_strat_3 mlmi_strat_4 mlmi_strat_5 mlmi_strat_6 mlmi_strat_7 mlmi_strat_8 log_w_bg log_hh log_sub mlmi_x17_age mlmi_x12_strat1 mlmi_x12_strat2 mlmi_x12_strat3"
    quietly _datasignature gl_subid gl_hhid gl_psu `ys' `xs'
    local input_signature "`r(datasignature)'"
    save "`base'"
    gen int mlmi_m=0
    save "`completed'"
    do code/stata/multilevel_mi.mata
    local perchain=`M'/`chains'
    forvalues c=1/`chains' {
        use "`base'", clear
        local seed=500000+`sim'+1000000*`c'
        set seed `seed'
        mata: MLMI = mlmi_initialize("`ys'","`xs'","gl_hhid","gl_psu", ///
            `ridge',(`hh_prior',`psu_prior',`person_prior'))
        if `c'==1 {
            local kept "`mlmi_kept_x'"
            preserve
            clear
            svmat double MLMI_scaling, names(col)
            gen str16 statistic=cond(_n==1,"observed_mean","observed_sd")
            export delimited using "`out'/sampler/scaling_`sim'.csv"
            restore
        }
        else assert "`kept'"=="`mlmi_kept_x'"
        mata: MLMI_FH=fopen("`out'/sampler/trace_`sim'_chain`c'.csv","w")
        mata: fput(MLMI_FH,mlmi_trace_header(MLMI))
        mata: mlmi_advance(MLMI,`burn',`c',`burn',`every',MLMI_FH)
        forvalues j=1/`perchain' {
            mata: mlmi_advance(MLMI,`between',`c',`burn',`every',MLMI_FH)
            use "`base'", clear
            sort gl_psu gl_hhid gl_subid
            mata: mlmi_store(MLMI)
            assert !missing(y_gfr_1,y_bmi_1,y_gfr_2,y_bmi_2,y_gfr_3,y_bmi_3)
            gen int mlmi_m=(`c'-1)*`perchain'+`j'
            append using "`completed'"
            save "`completed'", replace
            display "MLMI_DRAW sim=`sim' chain=`c' retained=`j'/`perchain'"
        }
        mata: fclose(MLMI_FH)
        mata: st_matrix("MLMI_acceptance",(MLMI.nc_accept,MLMI.hr_accept))
        preserve
        clear
        svmat double MLMI_acceptance, names(col)
        gen int chain=`c'
        export delimited using "`out'/sampler/acceptance_`sim'_chain`c'.csv"
        restore
    }
    use "`completed'", clear
    isid mlmi_m gl_subid
    * Import performs Stata's cross-imputation consistency checks.
    mi import flong, m(mlmi_m) id(gl_subid) ///
        imputed(y_bmi_2 y_gfr_2 y_bmi_3 y_gfr_3) clear
    mi convert wide, clear
    char _dta[mlmi_model] "joint6_psu_household_v1"
    char _dta[mlmi_status] "validation_is_in_sampler_json"
    char _dta[mlmi_input_signature] "`input_signature'"
    char _dta[mlmi_predictors] "`kept'"
    char _dta[mlmi_configuration] "M=`M';chains=`chains';burn=`burn';between=`between';every=`every';ridge=`ridge';priors=`hh_prior',`psu_prior',`person_prior'"
    quietly mi describe
    assert r(M)==`M'
    save "`candidate'"
    display "MLMI_GENERATION_COMPLETE sim=`sim' M=`M'; candidate only, no analysis fitted."
end

capture noisily _mlmi_generate `sim_id' `label' `M' `chains' `burn' `between' `every' `ridge' `hh_prior' `psu_prior' `person_prior'
local rc=_rc
capture log close mlmi
capture program drop _mlmi_generate
if `rc' exit `rc'
