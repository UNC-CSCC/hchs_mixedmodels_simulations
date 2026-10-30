* Verify the saved MI data against source observations; no models are fitted.
* Usage from diagnostic root: do code/stata/verify_mi_inputs.do <first> <last> <label> <M>
version 19
set more off
args first last label M
do code/stata/mi_rename_compat.do
log using "logs/oct2026/`label'/verify_mi_inputs_`first'_`last'.log", text replace name(micheck)
forvalues sim=`first'/`last' {
    import delimited "../../data/Oct2026/sample_miss00pct/samplemiss_`sim'.csv", clear varnames(1) case(lower)
    keep bgid hhid subid v_num y_gfr y_bmi miss_ind_mar
    rename y_gfr source_gfr
    rename y_bmi source_bmi
    rename miss_ind_mar source_mar
    tempfile source
    save "`source'"
    forvalues m=0/`M' {
        use "data/estimates/oct2026/`label'/imputations/mi_`sim'.dta", clear
        mi extract `m', clear
        reshape long y_gfr_ y_bmi_ x6_ mar_, i(gl_subid) j(v_num)
        merge 1:1 bgid hhid subid v_num using "`source'", assert(match) nogen
        assert mar_==source_mar
        assert y_gfr_==source_gfr & y_bmi_==source_bmi if source_mar==0
        if `m'==0 {
            assert missing(y_gfr_) & missing(y_bmi_) if source_mar==1
        }
        else {
            assert !missing(y_gfr_,y_bmi_)
        }
    }
    display "MI_INPUT_PASS sim=`sim': observed values preserved, flagged values masked in m=0, all M completed."
}
display "MI_INPUT_VERIFICATION_COMPLETE first=`first' last=`last' M=`M'"
log close micheck
