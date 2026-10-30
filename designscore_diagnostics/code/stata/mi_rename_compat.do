* Some local installations lack rename.ado. Use the installed primitive only
* when rename is unavailable. This fallback supports old_name new_name only.
* It exists only in this batch Stata session; no installation files are changed.
capture which rename
if !_rc exit
program define rename
    _rename `0'
end
display "MI compatibility: session-only rename delegates to built-in _rename."
