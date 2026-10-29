version 19
set more off
display c(stata_version)
sysdir
adopath
capture noisily which rename
capture noisily which meglm
capture noisily which mixed
capture noisily which egen
capture noisily which merge
exit, clear
