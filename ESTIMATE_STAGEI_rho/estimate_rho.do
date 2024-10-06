
**this code estimates STAGE I rho

clear
set more off, perm
set memory 100m
set obs 900
capture log close

* read raw data - change folder as needed
insheet using "//Users/kurtzv02/Library/CloudStorage/Dropbox/DN - uniform,pareto,2,6/data analysis codes/export/RawData/STAGE I/lottery_first_stage_2022-08-28.csv" , comma clear 

gen EV_lottery=0.5*playerx1_val + 0.5*playerx2_val

keep if sessioncode=="fxocegas"

rename playerslider_val ce
rename playerx1_val x1
rename playerx2_val x2


**mapping from stage 1 to stage 2, based on stage 1 ID and stage 1 participant code
gen id = 1 if participantcode=="brn25rv8"
replace id = 2 if participantcode=="11f3bm51"
replace id = 3 if participantcode=="p39qdtt5"
replace id = 4 if participantcode=="el6oyyvb"
replace id = 5 if participantcode=="wc3kni23"
replace id = 6 if participantcode=="3gzpsntj"
replace id = 7 if participantcode=="ewx2ztln"

gen session=16
gen uid = session*100+id

gen stage1_id = 1 if participantcode=="brn25rv8"
replace stage1_id = 2 if participantcode=="11f3bm51"
replace stage1_id = 3 if participantcode=="p39qdtt5"
replace stage1_id = 4 if participantcode=="el6oyyvb"
replace stage1_id = 5 if participantcode=="wc3kni23"
replace stage1_id = 6 if participantcode=="3gzpsntj"
replace stage1_id = 7 if participantcode=="ewx2ztln"


gen rho = .

**drop players who did not show up
drop if ce==.


**generate alpha
foreach i in 2 3 5 6  {
	disp "subject " `i'
nl (ce = (0.5*x1^{alpha=1.4} + 0.5*x2^{alpha=1.4})^(1/{alpha=1.4})) if id==`i'
scalar alpha_`i' = e(b)[1,1]
replace rho = alpha_`i' if id==`i'
}


preserve

bysort id: keep if _n==1

list id participantcode rho if rho!=. 

restore







