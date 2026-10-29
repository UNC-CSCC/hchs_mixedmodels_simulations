version 19
set more off
do scripts/stata/run_designscore_local.do 201 201 exact "data/design_frames/oct2026_subject_counts_201_smoke.csv"
display "DESIGNSCORE_201_SMOKE_OK"
exit, clear
