version 19
set more off
forvalues i=1/2 {
    mata: st_local("ds_loaded", strofreal(findexternal("ds_run()") != NULL))
    display "DESIGNSCORE_LOAD_CHECK iteration=`i' loaded=`ds_loaded'"
    if !`ds_loaded' quietly do scripts/stata/variance_designscore.mata
}
assert `ds_loaded' == 1
display "DESIGNSCORE_RELOAD_CHECK_OK"
exit, clear
