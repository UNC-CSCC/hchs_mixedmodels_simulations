version 19
set more off
args diagnostic_root
if `"`diagnostic_root'"' == "" exit 198
cd `"`diagnostic_root'"'
* Model-free check: only compile the Mata definitions.
quietly do code/stata/variance_designscore.mata
mata: st_local("loaded", strofreal(findexternal("ds_run()") != NULL))
assert `loaded' == 1
* The saved output MUST exist, so the model script must stop before any fit.
confirm file data/estimates/oct2026/exact/scen3_designscore_coef_1.csv
capture noisily do code/stata/scen3_designscore.do 1 oct2026 exact
local observed_rc = _rc
assert `observed_rc' == 602
display "RELOCATION_SMOKE_OK: Mata compiled, paths resolved, existing output protected; no models fitted."
exit, clear
