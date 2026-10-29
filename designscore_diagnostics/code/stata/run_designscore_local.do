version 19
set more off
* Invoke from designscore_diagnostics: do code/stata/run_designscore_local.do 1 200 exact
args first last mode counts_file
if "`mode'" == "" local mode exact
local first = real("`first'")
local last = real("`last'")
if missing(`first', `last') | `first'<1 | `last'<`first' exit 198
if !inlist("`mode'", "exact", "bernoulli") exit 198
local output "data/estimates/oct2026/`mode'"
capture mkdir data/estimates
capture mkdir data/estimates/oct2026
capture mkdir "`output'"
capture mkdir logs
capture mkdir logs/local_batches
tempname handle
local manifest "logs/local_batches/local_run_`first'_`last'.csv"
file open `handle' using "`manifest'", write text replace
file write `handle' "sim_id,return_code,status,elapsed_seconds" _n
file close `handle'
forvalues rep = `first'/`last' {
    local coef "`output'/scen3_designscore_coef_`rep'.csv"
    capture confirm file "`coef'"
    if !_rc {
        local rc 0
        local state "existing"
        local seconds 0
    }
    else {
        timer clear 99
        timer on 99
        capture noisily do code/stata/scen3_designscore.do `rep' oct2026 `mode' "" "" "" `"`counts_file'"'
        local rc = _rc
        timer off 99
        quietly timer list 99
        local seconds = r(t99)
        local state = cond(`rc'==0, "completed", "failed")
    }
    file open `handle' using "`manifest'", write text append
    file write `handle' "`rep',`rc',`state',`seconds'" _n
    file close `handle'
    display "LOCAL_DESIGNSCORE replicate=`rep' status=`state' rc=`rc' seconds=`seconds'"
    * Record an isolated nonconvergence and continue; stop on other errors.
    if `rc' != 0 & `rc' != 430 exit `rc'
}
display "LOCAL_DESIGNSCORE_BATCH_COMPLETE `first' `last' `mode'"
