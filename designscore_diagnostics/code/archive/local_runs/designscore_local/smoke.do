version 19
set more off
timer on 98
do scripts/stata/scen3_designscore.do 1 oct2026 exact "" "" "" "data/design_frames/oct2026_subject_counts_smoke.csv"
timer off 98
timer list 98
display "DESIGNSCORE_SMOKE_OK"
exit, clear
