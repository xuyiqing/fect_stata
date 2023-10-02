cap program drop fect
cap program drop cross_validation
cap program drop fect_counter
cap program drop placebo_Test
cap program drop carryover_Test
cap program drop permutation
cap program drop _gwtmean

program define fect,eclass
version 15.1

syntax varlist(min=1 max=1) [if], Treat(varlist min=1 max=1) ///
								  Unit(varlist min=1 max=1) ///
								  Time(varlist min=1 max=1) ///
							[ cov(varlist) ///
							weight(varlist max=1) ///
							  force(string) ///
							  degree(integer 3) ///
							  nknots(integer 3) ///
							  r(integer 4) ///
							  lambda(numlist) ///
							  nlambda(integer 10) ///
							  cv ///
							  kfolds(integer 10)  ///
							  cvtreat ///
							  cvnobs(integer 3) ///
							  method(string) ///
							  exit ///
							  alpha(real 0.05) ///
							  se ///
							  vartype(string) ///
							  nboots(integer 200) ///
							  tol(real 1e-5) ///
							  maxiterations(integer 2000) ///
							  seed(integer 100) ///
							  minT0(integer 5) ///
							  maxmissing(integer 0) ///
							  preperiod(integer -999999) ///
							  offperiod(integer 999999) ///
							  proportion(real 0.3) ///
							  wald ///
							  placeboTest ///
							  placeboperiod(integer 3) ///
							  carryoverTest ///
							  carryoverperiod(integer 3) ///
							  equiTest ///
							  permutation ///
							  title(string) ///
							  ylabel(string) ///
							  xlabel(string) ///
							  saving(string) ///
							  counterfactual ///
							]

set trace off

/*reghdfe*/
cap which reghdfe.ado
if _rc {
	di as error "reghdfe.ado required: {stata ssc install reghdfe,replace}"
	exit 111
}

cap which ftools.ado
if _rc {
	di as error "ftools.ado required: {stata ssc install ftools,replace}"
	exit 111
}

cap which _gwtmean.ado
if _rc {
	di as error "_gwtmean.ado required: {stata ssc install _gwtmean,replace}"
	exit 111
}

/* tokenize `varlist'
local y `1'
local treat `2'
macro shift
macro shift
local cov `*'
local varlist `y'

local unit `I'
local time `T' */

/* treat */
qui sum `treat',meanonly
if(r(min)!=0 | r(max)!=1){
	dis as err "treat() invalid; should only contain 0 and 1."
	exit
}

/* force */
local force = cond("`force'"=="","two-way","`force'")
if ("`force'"!="none" & "`force'"!="unit" & "`force'"!="time" & "`force'"!="two-way") {
	dis as err "force() invalid; choose one from none, unit, time and two-way."
	exit
	}
	
/* degree */
if(`degree'<1) {
	dis as err "degree() must be an integer greater than 1"
    exit
}

/* nknots */
if(`nknots'<3 | `nknots'>7) {
	dis as err "nknots() must be an integer between 3 and 7"
    exit
}

/* r */
if(`r'<1) {
	dis as err "r() must be an integer greater than 1"
    exit
}

/* lambda */
if("`lambda'"!="") {
	local smp = 0
	foreach lambda_grid in `lambda' {
		if(`lambda_grid'<=0){
				dis as err "lambda() invalid; elements in lambda must be positive."
				exit
		}
		local smp=`smp'+1
	}
	if(`smp'==1 & "`cv'"!=""){
		//dis as err "lambda() should have at least 2 values when cv is on."
		//exit
	}
}


/* nlambda */
if(`nlambda'<1) {
	dis as err "nlambda() must be an integer greater than 1"
    exit
}

/* cv */
if("`cv'"==""){
	local r = `r'
}
else{ 
	local nlambda = `nlambda'
}

/* kfold */
if(`kfolds'<3) {
	dis as err "k() must be an integer greater than 3"
    exit
}
local kfold=`kfolds'

/* cvnobs */							
if(`cvnobs'<1) {
	dis as err "cvnobs() must be an integer greater than 1"
    exit
}

/* method */
local method = cond("`method'"=="","fe","`method'")
if ("`method'"!="fe" & "`method'"!="ife" & "`method'"!="mc" & "`method'"!="bspline" & "`method'"!="polynomial" & "`method'"!="both") {
	dis as err "method() invalid; choose one from fe, ife, mc, bspline, polynomial and both."
	exit
	}

/* alpha */
if(`alpha'>0.5 | `alpha'<0) {
	dis as err "alpha() must be between 0 and 0.5"
    exit
}

/* vartype */
local vartype = cond("`vartype'"=="","bootstrap","`vartype'")
if ("`vartype'"!="bootstrap" & "`vartype'"!="jackknife") {
	dis as err "vartype() invalid; choose one from bootstrap and jackknife."
	exit
	}

/* nboots */
if(`nboots'<10) {
	dis as err "nboots() must be an integer greater than 10."
    exit
}	

/* tol */
if(`tol'<0) {
	dis as err "tol() must be larger than 0"
    exit
}

/* minT0 */
if(`minT0'<1) {
	dis as err "minT0() must be an integer greater than 1"
    exit
}

/* maxmissing */
if(`maxmissing'>0) {
	dis as res "Drop units whose number of missing values is larger than `maxmissing'."
}
if(`maxmissing'<0){
	dis as err "maxmissing() must be an integer greater than 0"
    exit
}

/* preperiod */
if(`preperiod'>0){
	dis as err "preperiod() must be an integer no larger than 0"
    exit
}

/* offperiod */
if(`offperiod'<0){
	dis as err "offperiod() must be an integer no less than 0"
    exit
}

/* Proportion */
if(`proportion'>=1 | `proportion'<0){
	dis as err "proportion() must be larger than 0 and less than 1"
    exit
}

/* placeboperiod */
if(`placeboperiod'<1) {
	dis as err "placeboperiod() must be an integer not less than 1"
    exit
}

/* carryoverperiod*/
if(`carryoverperiod'<1) {
    dis as err "carryoverperiod() must be an integer not less than 1"
    exit
}

/* saving */
if (strpos("`saving'",".") == 0) {
    local saving = "`saving'.png"
}  

/* EquiTest/PlaceboTest/CarryoverTest */
if("`se'"==""){
	local equiTest=""
	local placeboTest=""
	local carryoverTest=""
}

/*EquiTest*/
if("`equiTest'"!=""){
	local offperiod = 0
	dis "offperiod is set to 0 for Equivalence Test"
	if(`preperiod' == -999999){
	local preperiod = -2
	dis "Please state the number of preperiods for Equivalence Test. Default number is -2. "
	}

}


************************************Input Check End

	
/* save raw data */
tempfile raw
qui save `raw',replace

tempvar touse 
mark `touse' `if'
qui drop if `touse'==0
qui keep `varlist' `treat' `unit' `time' `cov' `weight' `touse'
qui drop if `unit' == .
qui drop if `time' == .

tempvar newid newtime
qui  egen `newid'=group(`unit')
sort `newid' `time'
qui  egen `newtime'=group(`time')

/* weight */
if("`weight'"!=""){
	tempvar sd_weight
	qui bys `newid': egen `sd_weight' = sd(`weight')
	qui sum `sd_weight',meanonly
	//if(r(mean)!=0){
	//	di as err "weights should be the same across periods for each unit."
	//	exit
	//}
}
else{
	qui cap drop `weight'
	tempvar weight
	qui gen `weight' = 1
}
	   
/* Check if the dataset is strongly balanced */	
qui tsset `newid' `newtime'					   
if "`r(balanced)'"!="strongly balanced"{
	di as res "Unbalanced Panel Data, fill the gap."
	tsfill,full
}
else{
	di as res "Balanced Panel Data"
}							   
sort `newid' `newtime'
qui mata: treat_fill("`treat'","`newid'","`newtime'") //fill the missing values in Treat
tempvar id
/* Check if the dataset is strongly balanced END */

/* Threshold of touse */
tempvar ftreat notreatnum notmissingy missingy
gen `ftreat'=0
qui replace `ftreat'=1 if `treat'==0 & `varlist'!=. //the number of training data
bys `newid': egen `notreatnum'=sum(`ftreat')
qui replace `touse'=0 if `notreatnum'<`minT0'

if(`maxmissing'>0){
	bys `newid': egen `notmissingy'=count(`varlist')
	qui sum `newtime',meanonly
	local TT = r(max)-r(min)
	qui gen `missingy'=`TT'-`notmissingy'
	qui replace `touse'=0 if `missingy'>`maxmissing'
}
qui drop if `touse'==0

if(r(N_drop)>0){
	di as res "Some treated units has too few pre-treatment periods; they are removed automatically."
}
/* Threshold of touse END*/

/* Re-Index */
tempvar newid2 newtime2
qui  egen `newid2'=group(`newid')
qui sort `newid2' `time'
qui  egen `newtime2'=group(`newtime')
qui sort `newid2' `newtime2'
qui drop `newid'
qui drop `newtime'
qui rename `newid2' `newid'
qui rename 	`newtime2' `newtime'
qui tsset `newid' `newtime'


/* period index: enter and exit*/
tempvar t_on
qui mata: gen_s("`treat'","`newid'","`newtime'","`t_on'")

tempvar t_off treat_rev
qui gen `treat_rev' = 1 - `treat'
qui mata: gen_s("`treat_rev'","`newid'","`newtime'","`t_off'")

qui sum `t_off'
if(r(N)==0){
	local hasrev = 0
}
else{
	local hasrev = 1
}

if(`hasrev'==1){
	di as res "Treatment has reversals."
}


/* Cross-Validation */
local lambda_size = wordcount("`lambda'")
local do_cv = 0
if("`method'"=="both"){
	local do_cv = 1
}
if("`cv'"!="" & "`method'"!="fe"){
	local do_cv = 1
}
if((`lambda_size'>1 | `lambda_size'==0)  & "`method'"=="mc"){
	local do_cv = 1
}

if(`do_cv'==1){
	cross_validation `varlist',treat(`treat') unit(`newid') time(`newtime') method(`method') t_on(`t_on') t_off(`t_off') ///
	weight(`weight') force(`force') cov(`cov') `cvtreat' cvnobs(`cvnobs') ///
	kfold(`kfold') nknots(`nknots') degree(`degree') ///
	r(`r') tol(`tol') maxiterations(`maxiterations') nlambda(`nlambda') lambda(`lambda') seed(`seed')
	
	local method_index=e(method)
	if(`method_index'==0){
		local method="fe"
		local lambda=0
		local r=0
	}
	else if(`method_index'==1){
		local method="ife"
		local r=e(optimal_parameter)
		local lambda=0
	}
	else if(`method_index'==2){
		local method="mc"
		local lambda=e(optimal_parameter)
		local r=0
	}
	else if(`method_index'==3){
		local method="bspline"
		local nknots=e(optimal_parameter)
		local lambda=0
		local r=0
	}
	else if(`method_index'==4){
		local method="polynomial"
		local degree=e(optimal_parameter)
		local lambda=0
		local r=0
	}
	
	local optimal_para=e(optimal_parameter) //to save
	tempname CV
	matrix `CV'=e(CV) //to save
	local min_mspe=e(min_mspe) //to save
}

if("`method'"!="mc"){
	local lambda=0
}

/* Estimate ATT & ATTs */
preserve
tempvar counter TE S ATTs ATTs_off

qui fect_counter `counter' `TE' `ATTs' `ATTs_off',outcome(`varlist') ///
treat(`treat') unit(`newid') time(`newtime') t_on(`t_on') ///
t_off(`t_off') method("`method'") force("`force'") ///
cov(`cov') nknots(`nknots') degree(`degree') r(`r') tol(`tol') ///
maxiterations(`maxiterations') lambda(`lambda') 

local judge_converge=e(converge)
if `judge_converge'==0 {
	dis as err "`method' can't converge."
    exit
}

//coef
tempname coef_output
if("`cov'"!=""){
	mat `coef_output' = e(coef)
	local colnum=colsof(`coef_output')
}
local cons_output = e(cons)


if("`counterfactual'"!=""){
	tempfile counter_out

	cap gen time_relative_to_treatment=`t_on'
	if _rc!=0 {
		cap replace time_relative_to_treatment=`t_on'
	}
	cap gen counterfactual_of_response=`counter'
	if _rc!=0 {
		cap replace counterfactual_of_response=`counter'
	}
	qui save `counter_out',replace

	qui drop time_relative_to_treatment
	qui drop counterfactual_of_response
	sort `newid' `newtime'

}

qui sum `t_on',meanonly
local min_s=r(min)
local max_s=r(max)

if(`preperiod'<`min_s'){
	local preperiod=`min_s'
}
else{
	local preperiod=`preperiod'
}

if(`offperiod'>`max_s'){
	local offperiod=`max_s'
}
else{
	local offperiod=`offperiod'
}

sort `newid' `newtime'
qui sum `TE' [iweight = `weight'] if `treat'==1 & `touse'==1 ,meanonly


local ATT=r(mean) //to save
local N_alltreat=r(N) //to save

tempname sim_nose
tempfile results_nose
qui postfile `sim_nose' s atts s_N using `results_nose',replace
forvalue s=`min_s'/`max_s' {
	qui sum `ATTs' [iweight = `weight'] if `t_on'==`s' & `touse'==1,meanonly
	if r(N)==0{
		//di in red "No observations for s=`s'"
		local s_index=`s'-`min_s'
		local att`s_index' = .
		local N`s_index'=0
	}
	else{
		local s_index=`s'-`min_s'
		local att`s_index' = r(mean)
		local N`s_index'=r(N)
	}
	post `sim_nose' (`s') (`att`s_index'') (`N`s_index'')
}
postclose `sim_nose'

if(`hasrev'==1){
	qui sum `t_off',meanonly
	local min_s_off=r(min)
	local max_s_off=r(max)

	if(`preperiod'<`min_s_off'){
		local preperiod_off=`min_s_off'
	}
	else{
		local preperiod_off=`preperiod'
	}

	if(`offperiod'>`max_s_off'){
		local offperiod_off=`max_s_off'
	}
	else{
		local offperiod_off=`offperiod'
	}

	tempname sim_nose_off
	tempfile results_nose_off
	qui postfile `sim_nose_off' s atts s_N using `results_nose_off',replace
	forvalue s=`min_s_off'/`max_s_off' {
		qui sum `ATTs_off' [iweight = `weight'] if `t_off'==`s' & `touse'==1,meanonly
		if r(N)==0{
			//di in red "No observations for s=`s'"
			local s_index=`s'-`min_s_off'
			local att_off`s_index' = .
			local N_off`s_index'=0
		}
		else{
			local s_index=`s'-`min_s_off'
			local att_off`s_index' = r(mean)
			local N_off`s_index'=r(N)
		}
		post `sim_nose_off' (`s') (`att_off`s_index'') (`N_off`s_index'')
	}
	postclose `sim_nose_off'

}

restore

/* Permutation Test */
if("`permutation'"!=""){
	permutation `varlist', treat(`treat') unit(`unit') time(`time') ///
	t_on(`t_on') cov(`cov') weight(`weight') ///
	force(`force') degree(`degree') nknots(`nknots') ///
	r(`r') lambda(`lambda') method(`method') nboots(`nboots') ///
	tol(`tol') maxiterations(`maxiterations') 

	local permutation_pvalue=e(permutation_pvalue)
}



/* Bootstrap/Jackknife */
if("`se'"!="" & "`placeboTest'"=="" & "`carryoverTest'"==""){ 
	tempname sim sim2 sim2_off sim3
	tempfile results results2 results2_off results3
	qui postfile `sim' ATT using `results',replace
	qui postfile `sim2' s ATTs using `results2',replace
	qui postfile `sim2_off' s ATTs using `results2_off',replace
	qui postfile `sim3' coef index using `results3',replace
	
	if("`vartype'"=="bootstrap"){
		di as txt "{hline}"	
		di as txt "Bootstrapping..."
		forvalue i=1/`nboots'{
			preserve
			tempvar counter TE S ATTs ATTs_off
			bsample, cluster(`newid') idcluster(boot_id) 
			drop `newid'
			rename boot_id `newid'
			
			capture qui fect_counter `counter' `TE' `ATTs' `ATTs_off',outcome(`varlist') treat(`treat') unit(`newid') time(`newtime') ///
			t_on(`t_on') t_off(`t_off') ///
			method("`method'") force("`force'") cov(`cov') nknots(`nknots') degree(`degree') r(`r') tol(`tol') ///
			maxiterations(`maxiterations') lambda(`lambda')
			
			if _rc!=0 {
				restore
				continue
			}
			
			/* coefs */
			local cons_output_boot = e(cons)
			post `sim3' (`cons_output_boot') (0)
			tempname coef_output_boot
			if("`cov'"!=""){
				mat `coef_output_boot' = e(coef)
				local colnum=colsof(`coef_output_boot')
				forval ii=1/`colnum'{
					post `sim3' (`coef_output_boot'[1,`ii']) (`ii')
				}
			}
			
			/* ATT */
			qui sum `TE' [iweight = `weight'] if `treat'==1 & `touse'==1,meanonly
			if r(mean)!=.{
				post `sim' (r(mean))
			}

			/* ATTs */
			forvalue si=`min_s'/`max_s' {
				qui sum `ATTs' [iweight = `weight'] if `t_on'==`si' & `touse'==1,meanonly
				if r(mean)!=.{
					post `sim2' (`si') (r(mean))
				}
			}

			/* ATTs with exit */
			if(`hasrev'==1){
				forvalue si=`min_s_off'/`max_s_off' {
					qui sum `ATTs_off' [iweight = `weight'] if `t_off'==`si' & `touse'==1,meanonly
					if r(mean)!=.{
						post `sim2_off' (`si') (r(mean))
					}
				}
			}
			if mod(`i',100)==0 {
				di as txt "ATT Estimation: Already Bootstrapped `i' Times"
			}
			restore
		}
		postclose `sim'
		postclose `sim2'   // all ATT_m calculate covariance F test. write meta. meta calculate covariance matrix
		postclose `sim2_off'
		postclose `sim3'
	}
	
	if("`vartype'"=="jackknife"){
		di as txt "{hline}"	
		di as txt "Jackknifing..."
		tempvar order order_mean order_rank
		bys `newid': gen `order'=runiform()
		bys `newid':  egen `order_mean'=mean(`order')
		sort `order_mean'
		egen `order_rank'=group(`order_mean')
		qui sum `newid',meanonly
		local max_id=r(max)
		forvalue i=1/`max_id'{
			preserve
			tempvar counter TE S ATTs ATTs_off
			qui drop if `order_rank'==`i'
			
			capture qui fect_counter `counter' `TE' `ATTs' `ATTs_off',outcome(`varlist') treat(`treat') unit(`newid') time(`newtime') ///
			t_on(`t_on') t_off(`t_off')  ///
			method("`method'") force("`force'") cov(`cov') nknots(`nknots') degree(`degree') r(`r') tol(`tol') ///
			maxiterations(`maxiterations') lambda(`lambda')
			
			if _rc!=0 {
				restore
				continue
			}

			/* coefs */
			local cons_output_jack = e(cons)
			post `sim3' (`cons_output_jack') (0)
			tempname coef_output_jack
			if("`cov'"!=""){
				mat `coef_output_jack' = e(coef)
				local colnum=colsof(`coef_output_jack')
				forval ii=1/`colnum'{
					post `sim3' (`coef_output_jack'[1,`ii']) (`ii')
				}
			}
			
			
			/* ATT */
			qui sum `TE' [iweight = `weight'] if `treat'==1 & `touse'==1,meanonly
			if r(mean)!=.{
				post `sim' (r(mean))
			}
		
			/* ATTs */
			forvalue si=`min_s'/`max_s' {
				qui sum `ATTs' [iweight = `weight'] if `t_on'==`si' & `touse'==1,meanonly
				if r(mean)!=.{
					post `sim2' (`si') (r(mean))
				}
			}
			if(`hasrev'==1){
				forvalue si=`min_s_off'/`max_s_off' {
					qui sum `ATTs' [iweight = `weight'] if `t_off'==`si' & `touse'==1,meanonly
					if r(mean)!=.{
						post `sim2_off' (`si') (r(mean))
					}
				}
			}
			//di as txt "Jackknife: Round `i' Finished"
			restore
		}
		postclose `sim'
		postclose `sim2'
		postclose `sim2_off'
		postclose `sim3'
	}
	
	
	//get estimation
	
	if("`vartype'"=="bootstrap"){
		//ATT
		preserve
		use `results',clear
		tempname ci ci_att
		qui sum ATT
		local ATTsd=r(sd) //to save
		local twosidel=100*`alpha'/2
		local twosideu=100*(1-`alpha'/2)
		local twoside=100-`alpha'*100
	
		qui _pctile ATT,p(`twosidel' ,`twosideu')
		local ATT_Lower_Bound=r(r1)
		local ATT_Upper_Bound=r(r2)
		matrix `ci_att'=(r(r1),r(r2)) //to save
		matrix rownames `ci_att' = "Confidence Interval" 
		matrix colnames `ci_att' = "Lower Bound" "Upper Bound"
	
		//ATTs
		tempname sim_plot
		tempfile results_plot
		qui postfile `sim_plot' s atts attsd att_lb att_ub s_N using `results_plot',replace
		qui use `results2',clear
		forvalue s=`min_s'/`max_s' {
			qui sum ATTs if s==`s'
			if r(N)<1 {
				continue
			}
			local attsd_temp=r(sd)
			qui _pctile ATTs if s==`s',p(`twosidel' ,`twosideu')
			matrix `ci'=(r(r1),r(r2))
			local s_index=`s'-`min_s'
			if("`vartype'"=="bootstrap"){
				post `sim_plot' (`s') (`att`s_index'') (`attsd_temp') (`ci'[1,1]) (`ci'[1,2]) (`N`s_index'')
			}
		}
		postclose `sim_plot'

		if(`hasrev'==1){
			//ATTs_off
			tempname sim_plot_off
			tempfile results_plot_off
			qui postfile `sim_plot_off' s atts attsd att_lb att_ub s_N using `results_plot_off',replace
			qui use `results2_off',clear
			forvalue s=`min_s_off'/`max_s_off' {
				qui sum ATTs if s==`s'
				if r(N)<1 {
					continue
				}
				local attsd_temp=r(sd)
				qui _pctile ATTs if s==`s',p(`twosidel' ,`twosideu')
				matrix `ci'=(r(r1),r(r2))
				local s_index=`s'-`min_s_off'
				if("`vartype'"=="bootstrap"){
					post `sim_plot_off' (`s') (`att_off`s_index'') (`attsd_temp') (`ci'[1,1]) (`ci'[1,2]) (`N_off`s_index'')
				}
			}
			postclose `sim_plot_off'			
		}

		//coef
		tempname sim_coef coef_ci
		tempfile results_coef
		qui postfile `sim_coef' coef sd p_value lower_bound upper_bound using `results_coef',replace
		qui use `results3',clear
		if("`cov'"==""){
			local colnum=0
		}
		forval ii=0/`colnum'{
			if(`ii'>0){
				local coef_tar=`coef_output'[1,`ii']
			}
			if(`ii'==0){
				local coef_tar=`cons_output'
			}
			qui sum coef if index==`ii'
			local coef_sd=r(sd)
			qui _pctile coef if index==`ii',p(`twosidel' ,`twosideu')
			matrix `coef_ci'=(r(r1),r(r2))
			local coef_p=round(2*(1-normal(abs(`coef_tar'/`coef_sd'))),0.001)
			
			post `sim_coef' (`coef_tar') (`coef_sd') (`coef_p') (`coef_ci'[1,1]) (`coef_ci'[1,2])
		}
		postclose `sim_coef'

		restore
	}
	
	
	if("`vartype'"=="jackknife"){
		//ATT
		preserve
		use `results',clear
		tempname ci ci_att
		qui sum ATT
		local ATTsd=r(sd)*sqrt(`max_id') //to save


		local tvalue = invttail(`max_id'-1,`alpha'/2)
		local ub_jack = `ATT'+`tvalue'*`ATTsd'
		local lb_jack = `ATT'-`tvalue'*`ATTsd'
		local ATT_Lower_Bound=`lb_jack'
		local ATT_Upper_Bound=`ub_jack'
		matrix `ci_att'=(`lb_jack',`ub_jack') //to save
		matrix rownames `ci_att' = "Confidence Interval" 
		matrix colnames `ci_att' = "Lower Bound" "Upper Bound"
	
		//ATTs
		tempname sim_plot
		tempfile results_plot
		qui postfile `sim_plot' s atts attsd att_lb att_ub s_N using `results_plot',replace
		qui use `results2',clear
		forvalue s=`min_s'/`max_s' {
			qui sum ATTs if s==`s'
			if r(N)<1 {
				continue
			}
			local attsd_temp=r(sd)*sqrt(`max_id')
			local s_index=`s'-`min_s'

			if("`vartype'"=="jackknife"){
				local tvalue = invttail(`max_id'-1,`alpha'/2)
				local ub_jack = `att`s_index''+`tvalue'*`attsd_temp'
				local lb_jack = `att`s_index''-`tvalue'*`attsd_temp'
				post `sim_plot' (`s') (`att`s_index'') (`attsd_temp') (`lb_jack') (`ub_jack') (`N`s_index'')
			}
		}
		postclose `sim_plot'

		if(`hasrev'==1){
			//ATTs_off
			tempname sim_plot_off
			tempfile results_plot_off
			qui postfile `sim_plot_off' s atts attsd att_lb att_ub s_N using `results_plot_off',replace
			qui use `results2_off',clear
			forvalue s=`min_s_off'/`max_s_off' {
				qui sum ATTs if s==`s'
				if r(N)<1 {
					continue
				}
				local attsd_temp=r(sd)*sqrt(`max_id')
				local s_index=`s'-`min_s_off'
				if("`vartype'"=="jackknife"){
					local tvalue = invttail(`max_id'-1,`alpha'/2)
					local ub_jack = `att_off`s_index''+`tvalue'*`attsd_temp'
					local lb_jack = `att_off`s_index''-`tvalue'*`attsd_temp'
					post `sim_plot_off' (`s') (`att_off`s_index'') (`attsd_temp') (`lb_jack') (`ub_jack') (`N_off`s_index'')
				}
			}
			postclose `sim_plot_off'			
		}

		//coef
		tempname sim_coef coef_ci
		tempfile results_coef
		qui postfile `sim_coef' coef sd p_value lower_bound upper_bound using `results_coef',replace
		qui use `results3',clear
		//list
		if("`cov'"==""){
			local colnum=0
		}
		forval ii=0/`colnum'{
			if(`ii'>0){
				local coef_tar=`coef_output'[1,`ii']
			}
			if(`ii'==0){
				local coef_tar=`cons_output'
			}
			qui sum coef if index==`ii'
			local coef_sd=r(sd)*sqrt(`max_id')
			local tvalue = invttail(`max_id'-1,`alpha'/2)
			local ub_jack = `coef_tar'+`tvalue'*`coef_sd'
			local lb_jack = `coef_tar'-`tvalue'*`coef_sd'

			matrix `coef_ci'=(`lb_jack',`ub_jack')
			local coef_p=round(2*(1-normal(abs(`coef_tar'/`coef_sd'))),0.001)
			post `sim_coef' (`coef_tar') (`coef_sd') (`coef_p') (`coef_ci'[1,1]) (`coef_ci'[1,2])
		}
		postclose `sim_coef'

		restore
	}
	
}


if("`se'"=="" & "`placeboTest'"=="" & "`carryoverTest'"=="" & "`equiTest'"==""){
	tempname sim_plot
	tempfile results_plot
	qui postfile `sim_plot' s atts s_N using `results_plot',replace
	
	forvalue s=`min_s'/`max_s' {
		local s_index=`s'-`min_s'
		post `sim_plot' (`s') (`att`s_index'') (`N`s_index'')
	}
	postclose `sim_plot'

	if(`hasrev'==1){
		tempname sim_plot_off
		tempfile results_plot_off
		qui postfile `sim_plot_off' s atts s_N using `results_plot_off',replace
		forvalue s=`min_s_off'/`max_s_off' {
			local s_index=`s'-`min_s_off'
			post `sim_plot_off' (`s') (`att_off`s_index'') (`N_off`s_index'')
		}
		postclose `sim_plot_off'
	}
}



/* Wald Test*/
if("`wald'"!=""){
	di as txt "{hline}"	
	di as txt "Wald Testing..."
	tempname sim_wald Fobs
	tempfile results_wald raw2
	qui postfile `sim_wald' SimuF using `results_wald',replace
	qui save `raw2',replace
	
	tempvar counter TE S ATTs ATTs_off
	qui fect_counter `counter' `TE' `ATTs' `ATTs_off',outcome(`varlist') ///
	treat(`treat') unit(`newid') time(`newtime') ///
	t_on(`t_on') t_off(`t_off')  ///
	method("`method'") force("`force'") cov(`cov') nknots(`nknots') degree(`degree') r(`r') tol(`tol') ///
	maxiterations(`maxiterations') lambda(`lambda')

	mata:waldF("`newid'", "`newtime'", "`TE'", "`ATTs'", "`t_on'", `preperiod', 0,"`touse'")
	local F_stat=r(F)
	
	forval i=1/`nboots' {
		preserve
		tempvar counter2 te2 s2 atts2 atts2_off w new_TE new_outcome
		qui gen `w'=runiformint(0,1)
		qui replace `w'=2*(`w'-0.5)
		qui gen `new_TE'=`w'*`TE'
		qui gen `new_outcome'=`counter'+`new_TE'
		
		qui fect_counter `counter2' `te2' `atts2' `atts2_off',outcome(`new_outcome') treat(`treat') unit(`newid') time(`newtime') ///
		t_on(`t_on') t_off(`t_off')  ///
		method("`method'") force("`force'") cov(`cov') nknots(`nknots') degree(`degree') r(`r') tol(`tol') ///
		maxiterations(`maxiterations') lambda(`lambda')
		
		sort `newid' `newtime'
		mata:waldF("`newid'", "`newtime'", "`te2'", "`atts2'", "`t_on'",  `preperiod', 0,"`touse'")
		
		post `sim_wald' (r(F))
		restore	
		if mod(`i',100)==0 {
			di as txt "Wald Testing: Already Simulated `i' Times"
		}	
	}
postclose `sim_wald'

qui use `results_wald',clear
qui sum SimuF,meanonly
qui count if SimuF>`F_stat'
local larger=r(N)
local WaldF=`F_stat'
local WaldP=`larger'/`nboots' //to save
local waldp_print=string(round(`WaldP',.001))
di as res "The p-value in wald test is `waldp_print'"

qui use `raw2',clear
}


/* Plot */

//draw plot0
if("`se'"=="" & "`placeboTest'"=="" & "`carryoverTest'"=="" & "`equiTest'"==""){
	preserve
	
	if("`exit'"!="" & `hasrev'==1){
		qui use `results_plot_off',clear
		local preperiod = `preperiod_off'
		local offperiod = `offperiod_off'
	}
	else{
		qui use `results_plot',clear
	}

	qui egen max_N=max(s_N)
	qui sum max_N,meanonly
	local max_N=r(mean)	
	local N_cutoff = `max_N'*`proportion'

	qui drop if s_N<`N_cutoff'
	qui drop if s<`preperiod'
	qui drop if s>`offperiod'

	//update preperiod and offperiod
	qui sum s, meanonly
	local preperiod = r(min)
	local offperiod = r(max)

	qui tsset s
	qui  egen max_atts=max(atts)

	local axis2_up=4*`max_N'
	qui sum max_atts,meanonly
	local max_atts=1.2*r(mean)
	if `max_atts'<0 {
		local max_atts=0
	}
	qui  egen min_atts=min(atts)
	qui sum min_atts,meanonly
	local min_atts=1.2*r(mean)
	if `min_atts'>0 {
		local min_atts=0
	}

	if(abs(`max_atts')>abs(`min_atts')){
		local min_atts=-`max_atts'*0.8
	}
	else{
		local max_atts=-`min_atts'*0.8
	}
	
	/* Title */
	if("`exit'"==""){
		local xlabel = cond("`xlabel'"=="","Time relative to the Treatment","`xlabel'")
	}
	else{
		local xlabel = cond("`xlabel'"=="","Time relative to the Exit of Treatment","`xlabel'")
	}
	local ylabel = cond("`ylabel'"=="","Average Treatment Effect","`ylabel'")
	local title = cond("`title'"=="","Estimated Average Treatment Effect","`title'") 
	if("`wald'"!=""){
		local note = "Wald p-value: `waldp_print'"
	}
	else{
		local note = " "
	}
	

	qui gen y0=0
	/*LSJ NEW*/
	twoway (scatter atts s, lcolor(black) yaxis(1) msize(2pt)) ///
	(bar s_N s,yaxis(2) barwidth(0.5) color(gray%50)), ///
	 xline(0, lcolor(gs6%50) lpattern(dash)) ///
	 yline(0, lcolor(gs6%50) lpattern(dash)) ///
	 yscale(axis(2) r(0 `axis2_up')) ///
	 yscale(axis(1) r(`min_atts' `max_atts')) ///
	 xscale(noextend) ///
	 xlabel(`preperiod' 0 (10) `offperiod', labsize(*0.75)) ///
	 xtick(`preperiod'(1)`offperiod') ///
	 ylabel(0 `max_N',axis(2) ) ///
	 ytitle(`ylabel',axis(1)) ///
	 ytitle("Num of observations",axis(2)) ///
	 xtitle(`xlabel') ///
	 title(`title') ///
	 text(`max_atts' `preperiod' "`note'",place(e)) ///
	 scheme(s2mono) ///
	 graphr(fcolor(white)) ///
	 plotregion(lcolor(black) lwidth(thin) margin(zero) ) ///
	 legend(order( 1 "ATT")) ///


	if("`saving'"!=""){
       cap graph export `saving',replace
	} 
	restore
}


//draw plot1
if("`se'"!="" & "`placeboTest'"=="" & "`carryoverTest'"=="" & "`equiTest'"==""){
	preserve
	if("`exit'"!="" & `hasrev'==1){
		qui use `results_plot_off',clear
		local preperiod = `preperiod_off'
		local offperiod = `offperiod_off'
	}
	else{
		qui use `results_plot',clear
	}

	qui egen max_N=max(s_N)
	qui sum max_N,meanonly
	local max_N=r(mean)	
	local N_cutoff = `max_N'*`proportion'

	qui drop if s_N<`N_cutoff'
	qui drop if s<`preperiod'
	qui drop if s>`offperiod'

	//update preperiod and offperiod
	qui sum s, meanonly
	local preperiod = r(min)
	local offperiod = r(max)

	qui tsset s

	local CIlevel=100*(1-`alpha')
	qui  egen max_atts=max(atts)
	local axis2_up=4*`max_N'
	qui sum max_atts,meanonly
	local max_atts=1.2*r(mean)
	if `max_atts'<0 {
		local max_atts=0
	}
	qui  egen min_atts=min(att_lb)
	qui sum min_atts,meanonly
	local min_atts=1.2*r(mean)
	if `min_atts'>0 {
		local min_atts=0
	}

	if(abs(`max_atts')>abs(`min_atts')){
		local min_atts=-`max_atts'*0.8
	}
	else{
		local max_atts=-`min_atts'*0.8
	}
	
	/* Title */
	if("`exit'"==""){
		local xlabel = cond("`xlabel'"=="","Time relative to the Treatment","`xlabel'")
	}
	else{
		local xlabel = cond("`xlabel'"=="","Time relative to the Exit of Treatment","`xlabel'")
	}	
	local ylabel = cond("`ylabel'"=="","Average Treatment Effect","`ylabel'")
	local title = cond("`title'"=="","Estimated Average Treatment Effect","`title'") 
	if("`wald'"!=""){
		local note = "Wald p-value: `waldp_print'"
	}
	else{
		local note = " "
	}

	 /*NEW changes LSJ */
	qui gen y0=0
	twoway (rcap att_lb att_ub s, color(black) lcolor(black) yaxis(1)) ///
	(scatter atts s, mcolor(black) yaxis(1) msize(2pt)) ///
	(bar s_N s,yaxis(2) barwidth(0.5) color(gray%50)), ///
	 xline(0.5, lcolor(gs6%50) lpattern(dash)) ///
	 yline(0, lcolor(gs6%50) lpattern(dash)) ///
	 yscale(axis(2) r(0 `axis2_up')) ///
	 yscale(axis(1) r(`min_atts' `max_atts')) ///
	 xscale(noextend) ///
	 xlabel(`preperiod' 0 (10) `offperiod', labsize(*0.75)) ///
	 xtick(`preperiod'(1)`offperiod') ///
	 ylabel(0 `max_N',axis(2) ) ///
	 ytitle(`ylabel',axis(1)) ///
	 ytitle("Num of observations",axis(2)) ///
	 xtitle(`xlabel') ///
	 title(`title') ///
	 text(`max_atts' `preperiod' "`note'",place(e)) ///
	 scheme(s2mono) ///
	 graphr(fcolor(white)) ///
	 plotregion(lcolor(black) lwidth(thin) margin(zero) ) ///
	 legend(order( 1 "ATT `CIlevel'% CI" 2 "ATT")) 
	 
	if("`saving'"!=""){
       cap graph export `saving',replace
	} 
	restore
}




/* Equivalence Test */
//draw plot2
// LSJ: Added F test
if("`se'"!="" & "`placeboTest'"=="" & "`carryoverTest'"=="" & "`equiTest'"!=""){
	di as txt "{hline}"	
	di as txt "Equivalence Test"
	qui reg `varlist' `cov' i.`newid' i.`newtime' [iweight = `weight'] if `treat'==0 & `touse'==1							   
	local ebound= 0.36*e(rmse)
	
	/* Etest is an one-side test*/
	local twosidel=100*`alpha'
	local twosideu=100*(1-`alpha')
	local twoside=100-`alpha'*100
	
	preserve
	tempname sim_etest
	tempfile results_etest
	qui postfile `sim_etest' s atts att_lb att_ub s_N using `results_etest',replace
	
	qui use `results2',clear
	local et_flag=0
	
	forvalue s=`min_s'/`max_s' {
		qui sum ATTs if s==`s',meanonly
		if r(N)<1 {
			continue
		}
	
		qui _pctile ATTs if s==`s',p(`twosidel' ,`twosideu')
		matrix `ci'=(r(r1),r(r2))
		local s_index=`s'-`min_s'
		post `sim_etest' (`s') (`att`s_index'') (`ci'[1,1]) (`ci'[1,2]) (`N`s_index'')
		
		if (`ci'[1,1]<(-`ebound') | `ci'[1,2]>(`ebound'))&(`s'<=0 & `s'>=`preperiod') { 
			di as res "Equivalence Test...Fail at s=`s'"
			local et_flag=1
		}
	}
	if(`et_flag'==1){
		di as res "Equivalence Test...Fail"
	}
	else{
		di as res "Equivalence Test...Pass"
	}
	postclose `sim_etest'
	
	/*F test*/
	local ft_flag=0
	qui use `results_etest', clear // still use the same one generated for TOST test.
	/*Get Ntr and m*/
	qui egen max_N=max(s_N)
	qui sum max_N,meanonly
	local max_N=r(mean)	
	local m_f = abs(`preperiod') 

	/*calculate F stats*/
	qui use `results2',clear
	qui mkmat ATTs if s==`preperiod', matrix(delta_m)
	local new_s_start = `preperiod'+1
	forvalue s=`new_s_start'/0{
		qui mkmat ATTs if s==`s', matrix(temp)
		qui matrix delta_m = delta_m, temp
	}
	/*Covariance matrix*/
	qui mata: delta_m = st_matrix("delta_m")
	qui	mata: vce_m = st_matrix("vce_m", variance(delta_m))

	qui mat list delta_m
	qui mat list vce_m // display the matrix
	/*Get F value*/
	qui mata: delta_vec = st_matrix("delta_vec", mean(delta_m))
	matrix F_mat = delta_vec*vce_m*delta_vec'
	scalar F = (`max_N'*(`max_N'-`m_f'-1)/((`max_N'-1)*(`m_f'+1))) * F_mat[1,1]

	di "F = " F
	local F_val = F
	local F_crit = invnFtail(`m_f',`max_N'-`m_f'-1,`m_f'*0.6,`alpha')
	if(`F_val' < `F_crit'){
		local ft_flag=1
	}
	if(`ft_flag'==1){
		di as res "Equivalence F-Test...Pass"
	}
	else{
		di as res "Equivalence F-Test...Fail"
	}	

	local FvalPrint = string(`F_val')
	/*F test ends*/

	/* Draw Plot */
	qui use `results_etest',clear

	qui egen max_N=max(s_N)
	qui sum max_N,meanonly
	local max_N=r(mean)	
	local N_cutoff = `max_N'*`proportion'
	qui drop if s_N<`N_cutoff'
	qui drop if s<`preperiod'
	qui drop if s>`offperiod'
	qui tsset s

	qui  egen max_atts=max(att_ub) if s<=0
	qui sum max_atts,meanonly
	qui gen yub=r(mean)
	local max_atts=2*r(mean)
	if `max_atts'<0 {
		local max_atts=0
	}
	if `max_atts'<`ebound' {
		local max_atts=`ebound'*2
	}
	qui  egen min_atts=min(att_lb) if s<=0
	qui sum min_atts,meanonly
	qui gen ylb=r(mean)
	local min_atts=2*r(mean)
	if `min_atts'>0 {
		local min_atts=0
	}
	if `min_atts'>-`ebound' {
		local min_atts=-`ebound'*2
	}
	
	local axis2_up=8*`max_N'
	local CIlevel=100*(1-`alpha')
	local Tlevel=100*(1-2*`alpha')
	qui gen y0=0
	qui gen y1=`ebound'
	qui gen y2=-`ebound'

	if(abs(`max_atts')>abs(`min_atts')){
		local min_atts=-`max_atts'*0.8
	}
	else{
		local max_atts=-`min_atts'*0.8
	}

	
	/* Title */
	local xlabel = cond("`xlabel'"=="","Time relative to the Treatment","`xlabel'")
	local ylabel = cond("`ylabel'"=="","Average Treatment Effect","`ylabel'")
	local title = cond("`title'"=="","Equivalence Test","`title'")
	
	/*NEW changes LSJ */
	twoway (rcap att_lb att_ub s, color(black)  lcolor(black) yaxis(1)) ///
	(scatter atts s, mcolor(black) yaxis(1) msize(2pt)) ///
	(line y1 y2 ylb yub s if s<=0,lcolor(blue blue green green) lpattern(dash dash dash dash) yaxis(1) lwidth(0.5 0.5 0.5 0.5)) ///
	(bar s_N s,yaxis(2) barwidth(0.5) color(gray%50)), ///
	 xline(0, lcolor(gs6%50) lpattern(dash)) ///
	 yline(0, lcolor(gs6%50) lpattern(dash)) ///
	 yscale(axis(2) r(0 `axis2_up')) ///
	 yscale(axis(1) r(`min_atts' `max_atts')) ///
	 xscale(noextend) ///
	 xlabel(`preperiod' 0 (10) `offperiod', labsize(*0.75)) ///
	 xtick(`preperiod'(1)`offperiod') ///
	 ylabel(0 `max_N',axis(2) ) ///
	 ytitle(`ylabel',axis(1)) ///
	 ytitle("Num of observations",axis(2)) ///
	 xtitle(`xlabel') ///
	 title(`title') ///
	 text(`max_atts' `preperiod' "`note'",place(e)) ///
	 text(`max_atts' `offperiod' "F-Test p-value: `FvalPrint'",place(sw)) ///
	 scheme(s2mono) ///
	 graphr(fcolor(white)) ///
	 plotregion(lcolor(black) lwidth(thin) margin(zero) ) ///
	legend(order(1 "ATT(`Tlevel'% CI)"  2 "ATT" 4 "Equiv.Bound" 6 "Min.Bound")) 
	
	
	if("`saving'"!=""){
       cap graph export `saving',replace
	} 
	restore

}

/* Placebo Test */
//draw plot3
if("`placeboTest'"!="" & "`se'"!="" & "`carryoverTest'"=="" & "`equiTest'"==""){

	placebo_Test `varlist',treat(`treat') unit(`newid') time(`newtime') ///
	cov(`cov') method("`method'") weight(`weight') proportion(`proportion') ///
	force("`force'") degree(`degree') nknots(`nknots') r(`r') lambda(`lambda') alpha(`alpha') vartype("`vartype'") ///
	nboots(`nboots') tol(`tol') maxiterations(`maxiterations') seed(`seed') minT0(`minT0') maxmissing(`maxmissing') ///
	preperiod(`preperiod') offperiod(`offperiod') placeboperiod(`placeboperiod') xlabel(`xlabel') ylabel(`ylabel') title(`title') saving(`saving')
	
	//to save
	local placebo_pvalue=e(placebo_pvalue)
	local placebo_ATT_mean=e(placebo_ATT_mean)
	local placebo_ATT_sd=e(placebo_ATT_sd)
	tempname placebo_CI
	matrix `placebo_CI' = e(placebo_CI)
	
}
 /*NEW LSJ*/
 /*Carryover Test*/
//draw plot4
if("`carryoverTest'"!="" & "`placeboTest'"=="" & "`se'"!="" & "`equiTest'"==""){
	//change preperiod and offperiod for exit
	/* if(`hasrev'==1){
		local preperiod = `preperiod_off'
		local offperiod = `offperiod_off'
	} */

	carryover_Test `varlist',treat(`treat') unit(`newid') time(`newtime') ///
	cov(`cov') method("`method'") weight(`weight') proportion(`proportion') ///
	force("`force'") degree(`degree') nknots(`nknots') r(`r') lambda(`lambda') alpha(`alpha') vartype("`vartype'") ///
	nboots(`nboots') tol(`tol') maxiterations(`maxiterations') seed(`seed') minT0(`minT0') maxmissing(`maxmissing') ///
	preperiod(`preperiod_off') offperiod(`offperiod_off') carryoverperiod(`carryoverperiod') xlabel(`xlabel') ylabel(`ylabel') title(`title') saving(`saving')
	
	//to save
	local carryover_pvalue=e(carryover_pvalue)
	local carryover_ATT_mean=e(carryover_ATT_mean)
	local carryover_ATT_sd=e(carryover_ATT_sd)
	tempname carryover_CI
	matrix `carryover_CI' = e(carryover_CI)
	
}

/*LSJ END*/

/* Storage */
ereturn clear
//CV
if( `do_cv'==1){
	ereturn scalar optimal_parameter=`optimal_para'
	ereturn scalar min_mspe=`min_mspe'
	ereturn matrix CV=`CV'
}
//ATT
//ereturn scalar ATT=`ATT'
//ereturn scalar n=`N_alltreat'

if("`se'"==""){
	preserve
	qui use `results_nose',clear
	qui rename atts ATT
	qui rename s_N n
	order s n ATT
	tempname ATTs_matrix
	qui mkmat s n ATT, matrix(`ATTs_matrix')
	ereturn matrix ATTs=`ATTs_matrix'
	if("`cov'"!=""){
		tempname coef_final_output
		matrix `coef_final_output'=(`cons_output')
		matrix `coef_final_output'=(`coef_final_output',`coef_output')
		matrix colnames `coef_final_output'=cons `cov'
		ereturn matrix coefs=`coef_final_output'
	}
	if("`cov'"==""){
		tempname coef_final_output
		matrix `coef_final_output'=(`cons_output')
		matrix colnames `coef_final_output'=cons
		ereturn matrix coefs=`coef_final_output'
	}

	tempname ATT_final_output
	matrix `ATT_final_output'=(`ATT',`N_alltreat')
	matrix colnames `ATT_final_output'=ATT N
	ereturn matrix ATT=`ATT_final_output'	

	restore
}

//ATTsd
if("`se'"!="" & "`placeboTest'"=="" & "`carryoverTest'"==""){
	//ereturn scalar ATTsd=`ATTsd'
	//ereturn matrix ATTCI=`ci_att'

	//ereturn scalar ATT_Lower_Bound=`ATT_Lower_Bound'
	//ereturn scalar ATT_Upper_Bound=`ATT_Upper_Bound'
	local ATT_pvalue=round(2*(1-normal(abs(`ATT'/`ATTsd'))),0.001)

	tempname ATT_final_output
	matrix `ATT_final_output'=(`ATT',`N_alltreat',`ATTsd',`ATT_Lower_Bound',`ATT_Upper_Bound',`ATT_pvalue')
	matrix colnames `ATT_final_output'=ATT N sd Lower_Bound Upper_Bound pvalue
	ereturn matrix ATT=`ATT_final_output'
}

//ATTs & coef
if("`se'"!="" & "`placeboTest'"=="" & "`carryoverTest'"==""){
	preserve
	qui use `results_plot',clear
	qui rename atts ATT
	qui rename attsd ATT_sd
	qui rename att_lb ATT_Lower_Bound
	qui rename att_ub ATT_Upper_Bound
	qui rename s_N n
	
	gen ATT_p_value=2*(1-normal(abs(ATT/ATT_sd)))
	order s n ATT ATT_sd ATT_p_value ATT_Lower_Bound ATT_Upper_Bound
	tempname ATTs_matrix
	qui mkmat s n ATT ATT_sd ATT_p_value ATT_Lower_Bound ATT_Upper_Bound, matrix(`ATTs_matrix')
	ereturn matrix ATTs=`ATTs_matrix'

	if(`hasrev'==1){
		qui use `results_plot_off',clear
		qui rename atts ATT
		qui rename attsd ATT_sd
		qui rename att_lb ATT_Lower_Bound
		qui rename att_ub ATT_Upper_Bound
		qui rename s_N n
		gen ATT_p_value=2*(1-normal(abs(ATT/ATT_sd)))
		order s n ATT ATT_sd ATT_p_value ATT_Lower_Bound ATT_Upper_Bound
		tempname ATTs_off_matrix
		qui mkmat s n ATT ATT_sd ATT_p_value ATT_Lower_Bound ATT_Upper_Bound, matrix(`ATTs_off_matrix')
		ereturn matrix ATTs_off=`ATTs_off_matrix'		
	}

	qui use `results_coef',clear
	tempname coefs_matrix_output
	qui mkmat coef sd p_value lower_bound upper_bound, matrix(`coefs_matrix_output')
	matrix rownames  `coefs_matrix_output'= cons `cov'
	ereturn matrix coefs=`coefs_matrix_output'
	restore
}
//permutation
if("`permutation'"!=""){
	ereturn scalar permutation_pvalue=`permutation_pvalue'
}

//Wald
if("`wald'"!=""){
	ereturn scalar wald_pvalue=`WaldP'
}

//Placebo
if("`placeboTest'"!=""){
	ereturn scalar placebo_pvalue=`placebo_pvalue'
	//ereturn scalar placebo_ATT=`placebo_ATT_mean'
	//ereturn scalar placebo_ATT_sd=`placebo_ATT_sd'
	//ereturn matrix placebo_ATT_CI=`placebo_CI'

	tempname placebo_ATT_final_output
	matrix `placebo_ATT_final_output'=(`placebo_ATT_mean',`placebo_ATT_sd',`placebo_CI'[1,1],`placebo_CI'[1,2],`placebo_pvalue')
	matrix colnames `placebo_ATT_final_output'=placebo_ATT sd Lower_Bound Upper_Bound pvalue
	ereturn matrix placebo_ATT=`placebo_ATT_final_output'

}

/*NEW LSJ*/
//Carryover
if("`carryoverTest'"!=""){
	ereturn scalar carryover_pvalue=`carryover_pvalue'

	tempname carryover_ATT_final_output
	matrix `carryover_ATT_final_output'=(`carryover_ATT_mean',`carryover_ATT_sd',`carryover_CI'[1,1],`carryover_CI'[1,2],`carryover_pvalue')
	matrix colnames `carryover_ATT_final_output'=carryover_ATT sd Lower_Bound Upper_Bound pvalue
	ereturn matrix carryover_ATT=`carryover_ATT_final_output'

}
/*LSJ END*/

//counterfactual
qui use `raw',clear	
if("`counterfactual'"!=""){
	//di as res "Counterfactual"
	qui use `counter_out',clear
	
}


		   
end

*********************************************************functions
program fect_counter,eclass
syntax newvarlist(min=4 max=4 gen) [if] , outcome(varlist min=1 max=1) ///
									treat(varlist min=1 max=1) ///
									unit(varlist min=1 max=1) ///
									time(varlist min=1 max=1) ///
									t_on(varlist min=1 max=1) ///
									t_off(varlist min=1 max=1) ///
									method(string) ///
									force(string) ///
									[ cov(varlist) ///
									  nknots(integer 3) ///
									  degree(integer 3) ///
									  r(integer 2) ///
									  tol(real 1e-5) ///
									  maxiterations(integer 5000) ///
									  lambda(real 0.5) ///
									  t_on_original(varlist max=1) ///
									  t_off_original(varlist max=1) ///
									]

set trace off
tempvar touse 
mark `touse' `if'

tempvar newid newtime
qui  egen `newid'=group(`unit')
sort `newid' `time'
qui  egen `newtime'=group(`time')

/* Normalization */
tempvar norm_y
//qui sum `outcome' if `treat'==0 & `touse'==1
//local mean_y=r(mean)
//local sd_y=r(sd)
//qui gen `norm_y'=(`outcome'-`mean_y')/`sd_y'
qui gen `norm_y'=`outcome'

/* Regression/Training */

if("`method'"=="fe"){
	tempvar FE1 FE2 yearFE unitFE counterfactual
	if("`force'"=="two-way"){
		qui reghdfe `norm_y' `cov' if `treat'==0 & `touse'==1,absorb(`FE1'=`newid' `FE2'=`newtime') ///
		  dof(none)   
	}
	if("`force'"=="unit"){
		qui reghdfe `norm_y' `cov' if `treat'==0 & `touse'==1,absorb(`FE1'=`newid')   ///
		dof(none)   
		qui gen `FE2'=0
	}
	if("`force'"=="time"){
		qui reghdfe `norm_y' `cov' if `treat'==0 & `touse'==1,absorb(`FE2'=`newtime')   ///
		dof(none)   
		qui gen `FE1'=0
	}
	if("`force'"=="none"){
		qui reghdfe `norm_y' `cov' if `treat'==0 & `touse'==1,noabsorb ///
		dof(none)   
		qui gen `FE2'=0
		qui gen `FE1'=0
	}

	qui bys `newtime': egen `yearFE'=mean(`FE2')
	qui bys `newid': egen `unitFE'=mean(`FE1')
	drop `FE1' `FE2'
	qui gen double `counterfactual'=_b[_cons]+`unitFE'+`yearFE' 
	if "`cov'"!="" {
		tempname coefficient
		matrix `coefficient'=e(b)
		local colnum=colsof(`coefficient')-1
		forval i=1/`colnum'{
			local p=word("`cov'",`i')
			qui replace `counterfactual'=`counterfactual'+`p'*`coefficient'[1,`i']
		}
	}
	local Converge=1
	tempname coef_output
	if( "`cov'"!=""){
		forval i=1/`colnum'{
			if(`i'==1){
				matrix `coef_output'=(`coefficient'[1,`i'])
			}
			else{
				matrix `coef_output'=(`coef_output',`coefficient'[1,`i'])
			}
		}
	}
	local mu_output=_b[_cons]
}

if("`method'"=="bspline"){
	tempvar trend
	tempvar FE1 FE2 yearFE unitFE counterfactual FE
	qui mkspline `trend' = `newtime',cubic nknots(`nknots')
	local run = `nknots'-1
	forvalue i=1/`run'{
		local trend_all `trend_all' `trend'`i'
		tempvar `FE'Slope`i'
	}
	
	if("`force'"=="two-way"){
		qui reghdfe `norm_y' `cov' if `treat'==0 & `touse'==1,absorb(`FE1'=`newid' `FE2'=`newtime' `FE'=i.`newid'#c.(`trend_all'))   ///
		dof(none)   
	}
	if("`force'"=="unit"){
		qui reghdfe `norm_y' `cov' if `treat'==0 & `touse'==1,absorb(`FE1'=`newid' `FE'=i.`newid'#c.(`trend_all'))   ///
		dof(none)   
		qui gen `FE2'=0
	}
	if("`force'"=="time"){
		qui reghdfe `norm_y' `cov' if `treat'==0 & `touse'==1,absorb(`FE2'=`newtime' `FE'=i.`newid'#c.(`trend_all'))   ///
		dof(none)   
		qui gen `FE1'=0
	}
	if("`force'"=="none"){
		qui reghdfe `norm_y' `cov' if `treat'==0 & `touse'==1,absorb(`FE'=i.`newid'#c.(`trend_all'))   ///
		dof(none)   
		qui gen `FE2'=0
		qui gen `FE1'=0
	}
	
	qui bys `newtime': egen `yearFE'=mean(`FE2')
	qui bys `newid': egen `unitFE'=mean(`FE1')
	drop `FE1' `FE2'
	
	if("`force'"!="none"){
		qui gen `counterfactual'=_b[_cons]+`unitFE'+`yearFE' 
	}
	else{
		qui gen `counterfactual'=0
	}
	
	
	if "`cov'"!="" {
		tempname coefficient
		matrix `coefficient'=e(b)
		if("`force'"!="none"){
			local colnum=colsof(`coefficient')-1
		}
		else{
			local colnum=colsof(`coefficient')
		}
		forval i=1/`colnum'{
			local p=word("`cov'",`i')
			qui replace `counterfactual'=`counterfactual'+`p'*`coefficient'[1,`i']
		}
	}
	
	forvalue i=1/`run'{
		qui bys `newid': egen `unitFE'`i'=mean(`FE'Slope`i')
		qui replace `counterfactual'=`counterfactual'+`unitFE'`i'*`trend'`i'
	}
	local Converge=1
	tempname coef_output
	if( "`cov'"!=""){
		forval i=1/`colnum'{
			if(`i'==1){
				matrix `coef_output'=(`coefficient'[1,`i'])
			}
			else{
				matrix `coef_output'=(`coef_output',`coefficient'[1,`i'])
			}
		}
	}
	local mu_output=_b[_cons]
}

if("`method'"=="polynomial"){
	tempvar trend
	tempvar FE1 FE2 yearFE unitFE counterfactual FE
	forvalue i=1/`degree'{
		qui gen `trend'`i' = `newtime'^`i'
		local trend_all `trend_all' `trend'`i'
	}
	
	
	if("`force'"=="two-way"){
		qui reghdfe `norm_y' `cov' if `treat'==0 & `touse'==1,absorb(`FE1'=`newid' `FE2'=`newtime' `FE'=i.`newid'#c.(`trend_all'))   ///
		dof(none)   
	}
	if("`force'"=="unit"){
		qui reghdfe `norm_y' `cov' if `treat'==0 & `touse'==1,absorb(`FE1'=`newid' `FE'=i.`newid'#c.(`trend_all'))   ///
		dof(none)   
		qui gen `FE2'=0
	}
	if("`force'"=="time"){
		qui reghdfe `norm_y' `cov' if `treat'==0 & `touse'==1,absorb(`FE2'=`newtime' `FE'=i.`newid'#c.(`trend_all'))   ///
		dof(none)   
		qui gen `FE1'=0
	}
	if("`force'"=="none"){
		qui reghdfe `norm_y' `cov' if `treat'==0 & `touse'==1,absorb(`FE'=i.`newid'#c.(`trend_all'))   ///
		dof(none)   
		qui gen `FE2'=0
		qui gen `FE1'=0
	}
	
	qui bys `newtime': egen `yearFE'=mean(`FE2')
	qui bys `newid': egen `unitFE'=mean(`FE1')
	drop `FE1' `FE2'
	
	if("`force'"!="none"){
		qui gen `counterfactual'=_b[_cons]+`unitFE'+`yearFE' 
	}
	else{
		qui gen `counterfactual'=0
	}
	
	if "`cov'"!="" {
		tempname coefficient
		matrix `coefficient'=e(b)
		if("`force'"!="none"){
			local colnum=colsof(`coefficient')-1
		}
		else{
			local colnum=colsof(`coefficient')
		}
		forval i=1/`colnum'{
			local p=word("`cov'",`i')
			qui replace `counterfactual'=`counterfactual'+`p'*`coefficient'[1,`i']
		}
	}
	
	forvalue i=1/`degree'{
		qui bys `newid': egen `unitFE'`i'=mean(`FE'Slope`i')
		qui replace `counterfactual'=`counterfactual'+`unitFE'`i'*`trend'`i'
	}
	local Converge=1
	tempname coef_output
	if( "`cov'"!=""){
		forval i=1/`colnum'{
			if(`i'==1){
				matrix `coef_output'=(`coefficient'[1,`i'])
			}
			else{
				matrix `coef_output'=(`coef_output',`coefficient'[1,`i'])
			}
		}
	}
	local mu_output=_b[_cons]
}

if("`method'"=="ife"){
	tempvar FE1 FE2 yearFE unitFE cons counterfactual p_treat
	if("`force'"=="two-way"){
		qui reghdfe `norm_y' `cov' if `treat'==0 & `touse'==1,absorb(`FE1'=`newid' `FE2'=`newtime')   ///
		dof(none)   
	}
	if("`force'"=="unit"){
		qui reghdfe `norm_y' `cov' if `treat'==0 & `touse'==1,absorb(`FE1'=`newid')   ///
		dof(none)   
		qui gen `FE2'=0
	}
	if("`force'"=="time"){
		qui reghdfe `norm_y' `cov' if `treat'==0 & `touse'==1,absorb(`FE2'=`newtime')   ///
		dof(none)   
		qui gen `FE1'=0
	}
	if("`force'"=="none"){
		qui reghdfe `norm_y' `cov' if `treat'==0 & `touse'==1,noabsorb ///
		dof(none)   
		qui gen `FE2'=0
		qui gen `FE1'=0
	}
	qui gen `cons'=_b[_cons]
	qui bys `newtime': egen `yearFE'=mean(`FE2')
	qui bys `newid': egen `unitFE'=mean(`FE1')
	drop `FE1' `FE2'
	
	sort `newid' `newtime'
	qui gen `p_treat'=0 
	qui replace `p_treat'=1 if `treat'==1 | `touse'==0 | `norm_y'==.
	qui mata: ife("`norm_y'","`cov'","`p_treat'","`newid'","`newtime'","`force'",`r',`tol',`maxiterations',"`counterfactual'","`unitFE'","`yearFE'","`cons'")
	if(r(stop)==1){
		local Converge=1
	}
	else{
		local Converge=0
		di as res "IFE can't converge."
	}
	tempname coef_output
	if "`cov'"!=""{
		matrix `coef_output'=r(coef)
	}
	local mu_output=r(cons)

}

if("`method'"=="mc"){
	tempvar y FE1 FE2 yearFE unitFE cons counterfactual p_treat
	qui gen `y'=`outcome'
	if("`force'"=="two-way"){
		qui reghdfe `y' `cov' if `treat'==0 & `touse'==1,absorb(`FE1'=`newid' `FE2'=`newtime')   ///
		dof(none)   
	}
	if("`force'"=="unit"){
		qui reghdfe `y' `cov' if `treat'==0 & `touse'==1,absorb(`FE1'=`newid')   ///
		dof(none)   
		qui gen `FE2'=0
	}
	if("`force'"=="time"){
		qui reghdfe `y' `cov' if `treat'==0 & `touse'==1,absorb(`FE2'=`newtime')   ///
		dof(none)   
		qui gen `FE1'=0
	}
	if("`force'"=="none"){
		qui reghdfe `y' `cov' if `treat'==0 & `touse'==1,noabsorb ///
		dof(none)   
		qui gen `FE2'=0
		qui gen `FE1'=0
	}
	qui gen `cons'=_b[_cons]
	qui bys `newtime': egen `yearFE'=mean(`FE2')
	qui bys `newid': egen `unitFE'=mean(`FE1')
	drop `FE1' `FE2'
	
	sort `newid' `newtime'
	qui gen `p_treat'=0 
	qui replace `p_treat'=1 if `treat'==1 | `touse'==0 | `y'==.
	qui mata: mc("`y'","`cov'","`p_treat'","`newid'","`newtime'","`force'",`lambda',`tol',`maxiterations',"`counterfactual'","`unitFE'","`yearFE'","`cons'")
	if(r(stop)==1){
		local Converge=1
	}
	else{
		local Converge=0
		di as res "MC can't converge."
	}
	
	tempname coef_output
	if "`cov'"!=""{
		matrix `coef_output'=r(coef)
	}
	local mu_output=r(cons)
}

if("`cov'"!=""){
	matrix colnames `coef_output' = `cov'
}

token `varlist'
qui replace `1' = `counterfactual'

tempvar treat_effect
qui gen double `treat_effect'=`outcome'-`counterfactual' 
qui replace `2' = `treat_effect'
sort `unit' `time'

tempvar treatsum targettime target_on ATTs target_off ATTs_off
//qui mata:gen_s("`treat'","`newid'","`newtime'","`target'")
gen `target_on'=`t_on'
if("`t_on_original'"!=""){
	replace `target_on'=`t_on_original'
}
qui bys `target_on': egen `ATTs'= mean(`treat_effect') if `target_on'!=. & `touse'==1
qui replace `3' = `ATTs' if `touse'==1

gen `target_off'=`t_off'
if("`t_off_original'"!=""){
	replace `target_off'=`t_off_original'
}
qui bys `target_off': egen `ATTs_off'= mean(`treat_effect') if `target_off'!=. & `touse'==1
qui replace `4' = `ATTs_off' if `touse'==1

qui sort `newid' `newtime'
ereturn scalar converge=`Converge'
ereturn scalar cons= `mu_output'
if("`cov'"!=""){
	ereturn matrix coef= `coef_output'
}

end

program cross_validation,eclass
syntax varlist(min=1 max=1), treat(varlist min=1 max=1) ///
						     unit(varlist min=1 max=1) ///
						     time(varlist min=1 max=1) ///
							 t_on(varlist min=1 max=1) ///
							 t_off(varlist min=1 max=1) ///
						     method(string) ///
							 force(string) ///
							[ cov(varlist) ///
							  weight(varlist max=1) ///
							  cvtreat ///
							  cvnobs(integer 3) ///
							  kfold(integer 10) ///
							  nknots(integer 3) ///
							  degree(integer 3) ///
							  r(integer 3) ///
							  tol(real 1e-4) ///
							  maxiterations(integer 5000) ///
							  nlambda(integer 10) ///
							  lambda(numlist ascending) ///
							  seed(integer 123) ///
						    ]

set trace off
di as txt "{hline}"	
di as txt "Cross Validation..."
tempvar newid newtime touse
qui gen `touse'=1
qui  egen `newid'=group(`unit')
qui  egen `newtime'=group(`time')
sort `newid' `newtime'

tempvar validation
tempname sim final_sim final_sim2
tempfile results final_results final_results2	
gen `validation'=0

/* Only for treated values */
if("`cvtreat'"!=""){
	tempvar iftreat obs_num
	tempfile tempsave
	qui gen `obs_num' = _n
	qui save `tempsave',replace
	qui bys `newid':  egen `iftreat'=mean(`treat')
	qui keep if `iftreat'>0
	mata:val_gen("`treat'", "`newid'", "`newtime'", "`touse'", ///
				 `kfold',`cvnobs',"`validation'",`seed')
	qui keep `obs_num' `validation'
	qui merge 1:1 `obs_num' using `tempsave', nogen nol nonote norep
	qui drop `obs_num'
	qui replace `validation' = 0 if `validation'==.
	sort `newid' `newtime'
}
else{
	mata: val_gen("`treat'", "`newid'", "`newtime'", "`touse'", ///
	`kfold',`cvnobs',"`validation'",`seed')
}

/*
if("`cvtreat'"!=""){ 
	tempvar iftreat
	qui bys `newid':  egen `iftreat'=mean(`treat')
	qui replace `validation'=0 if `iftreat'==0
}
*/

if("`method'"=="ife" | "`method'"=="both"){
	qui postfile `final_sim' r mspe using `final_results',replace
	forval fnum=0/`r' {
		//di as res "Factor: `fnum'"
		qui postfile `sim' N spe using `results',replace
		if `fnum'==0 {
			forval i=1/`kfold' {
				preserve
				tempvar counter spe TE S ATTs ATTs_off
				qui fect_counter `counter' `TE' `ATTs' `ATTs_off' ///
				if `validation'!=`i' & `touse'==1,outcome(`varlist') ///
				treat(`treat') unit(`newid') time(`newtime') t_on(`t_on') ///
				t_off(`t_off') cov(`cov') method("fe") force("`force'") 
				
				if(e(converge)==0){
					restore
					continue
				}

				qui gen `spe'=(`varlist'-`counter')^2 if `validation'==`i' & `touse'==1
				qui sum `spe' [iweight=`weight'] if `validation'==`i' & `touse'==1,meanonly
				if r(N)!=0 {
					post `sim' (r(N)) (r(mean))
				}
				restore
			}
		postclose `sim'

		preserve
		use `results',clear
		tempvar weight_mean
		egen `weight_mean'=wtmean(spe),weight(N)
		qui sum `weight_mean',meanonly
		post `final_sim' (0) (r(mean))
		local mspe = string(round(r(mean),.0001))
		di as txt "fe r=0 force=`force' mspe=`mspe'"
		restore
		continue
		}

		/* ife */
		forval i=1/`kfold' {
			preserve
			tempvar counter spe TE S ATTs ATTs_off

			qui fect_counter `counter' `TE' `ATTs' `ATTs_off' ///
			if `validation'!=`i' & `touse'==1,outcome(`varlist') ///
			treat(`treat') unit(`newid') time(`newtime') t_on(`t_on') ///
			t_off(`t_off') cov(`cov') method("ife") force("`force'") r(`fnum')  ///
			tol(`tol') maxiterations(`maxiterations') 
			
			if(e(converge)==0){
				restore
				continue
			}
				
			qui gen `spe'=(`varlist'-`counter')^2 if `validation'==`i' & `touse'==1
			qui sum `spe' [iweight=`weight'] if `validation'==`i' & `touse'==1,meanonly

			if r(N)!=0 {
				post `sim' (r(N)) (r(mean))
			}

			restore
			//di as res "Fold `i' finished"
		}
		postclose `sim'


		preserve
		use `results',clear
		tempvar weight_mean
		 egen `weight_mean'=wtmean(spe),weight(N)
		qui sum `weight_mean',meanonly
		local mspe = string(round(r(mean),.0001))
		post `final_sim' (`fnum') (r(mean))
		di as txt "ife r=`fnum' force=`force' mspe=`mspe'"
		restore
	}
	postclose `final_sim'
}

if("`method'"=="mc" | "`method'"=="both"){

	qui postfile `final_sim2' lambda lambda_norm mspe using `final_results2',replace
	//get lambda list
	tempvar FE1 FE2 e
	if("`force'"=="two-way"){
		qui reghdfe `varlist' `cov' if `treat'==0 & `touse'==1, ab(`FE1'=`newid' `FE2'=`newtime') resid ///
		dof(none)   
	}
	if("`force'"=="unit"){
		qui reghdfe `varlist' `cov' if `treat'==0 & `touse'==1, ab(`FE1'=`newid') resid ///
		dof(none)   
	}
	if("`force'"=="unit"){
		qui reghdfe `varlist' `cov' if `treat'==0 & `touse'==1, ab(`FE2'=`newtime') resid ///
		dof(none)   
	}
	if("`force'"=="none"){
		qui reghdfe `varlist' `cov' if `treat'==0 & `touse'==1, noabsorb resid ///
		dof(none)   
	}
		
	qui predict `e',residual
		
	qui mata: eig_v("`e'","`treat'","`newid'","`newtime'")
		
	local lambda_max=log10(r(max_sing))
	local lambda_by = 3/(`nlambda' - 2)
	local start_lambda= 10^(`lambda_max')
		
	if("`lambda'"==""){
		if(`nlambda'<3){
			di as err "nlambda should be larger than 2."
			exit
		}
	
		
		local lambda `start_lambda'
		
		forval num=2/`nlambda'{
			local lambda_add=10^(`lambda_max' - (`num' - 1)*`lambda_by')
			local lambda `lambda' `lambda_add'   
		}		
					
	}
		
	foreach lambda_grid in `lambda' {
		local lambda_grid_norm = `lambda_grid'/`start_lambda'
		qui postfile `sim' N spe using `results',replace

		forval i=1/`kfold' {
			preserve
			tempvar counter spe TE S ATTs ATTs_off
			
			qui fect_counter `counter' `TE' `ATTs' `ATTs_off' if `validation'!=`i' & `touse'==1,outcome(`varlist') ///
			treat(`treat') unit(`newid') time(`newtime') t_on(`t_on') ///
			t_off(`t_off') cov(`cov') method("mc") force(`force') ///
			tol(`tol') maxiterations(`maxiterations') lambda(`lambda_grid')
			
			if(e(converge)==0){
				restore
				continue
			}

			qui gen `spe'=(`varlist'-`counter')^2  if `validation'==`i' & `touse'==1
			qui sum `spe' [iweight=`weight'] if `validation'==`i' & `touse'==1 ,meanonly
			if r(N)!=0 {
				post `sim' (r(N)) (r(mean))
			}
				restore
		}
		postclose `sim'

		preserve
		use `results',clear
		tempvar weight_mean
		egen `weight_mean'=wtmean(spe),weight(N)
		qui sum `weight_mean' ,meanonly
		post `final_sim2' (`lambda_grid') (`lambda_grid_norm') (r(mean))
		local mspe = string(round(r(mean),.0001))
		local lambda_print = string(round(`lambda_grid',.0001))
		local lambda_norm_print = string(round(`lambda_grid'/`start_lambda',.0001))
		di as txt "mc: lambda=`lambda_print' lambda.norm=`lambda_norm_print' mspe=`mspe'"
		restore
	}
	postclose `final_sim2'
}

if("`method'"=="bspline"){
	qui postfile `final_sim' nknots mspe using `final_results',replace
	if(`nknots'>7 | `nknots'<3){
		dis as err "nkots must be between 3 and 7"
		exit
	}
	
	forval knot=3/`nknots' {
		qui postfile `sim' N spe using `results',replace

		forval i=1/`kfold' {
			preserve
			tempvar counter spe TE S ATTs ATTs_off
			qui fect_counter `counter' `TE'  `ATTs' `ATTs_off' if `validation'!=`i' & `touse'==1,outcome(`varlist') ///
			treat(`treat') unit(`newid') time(`newtime') t_on(`t_on') ///
			t_off(`t_off') cov(`cov') method("bspline") nknots(`knot') force(`force')

			qui gen `spe'=(`varlist'-`counter')^2 if `validation'==`i' & `touse'==1
			qui sum `spe' [iweight=`weight'] if `validation'==`i' & `touse'==1,meanonly
			if r(N)!=0 {
				post `sim' (r(N)) (r(mean))
			}
				restore
		}
		postclose `sim'

		preserve
		use `results',clear
		tempvar weight_mean
		 egen `weight_mean'=wtmean(spe),weight(N)
		qui sum `weight_mean',meanonly
		post `final_sim' (`knot') (r(mean))
		local mspe = string(round(r(mean),.0001))
		di as txt "bspline nknots=`knot' mspe=`mspe'"
		restore
	}
	
	
	postclose `final_sim'
}

if("`method'"=="polynomial"){
	qui postfile `final_sim' degree mspe using `final_results',replace
	if(`degree'<1){
		dis as err "degree must be a positive integer."
		exit
	}
	
	forval dg=1/`degree' {
		qui postfile `sim' N spe using `results',replace

		forval i=1/`kfold' {
			preserve
			tempvar counter spe TE S ATTs ATTs_off
			qui fect_counter `counter' `TE' `ATTs' `ATTs_off' if `validation'!=`i' & `touse'==1,outcome(`varlist') ///
			treat(`treat') unit(`newid') time(`newtime') t_on(`t_on') ///
			t_off(`t_off') cov(`cov') method("polynomial") degree(`dg') force(`force')

			qui gen `spe'=(`varlist'-`counter')^2 if `validation'==`i' & `touse'==1
			qui sum `spe' [iweight=`weight'] if `validation'==`i' & `touse'==1,meanonly
			if r(N)!=0 {
				post `sim' (r(N)) (r(mean))
			}
				restore
		}
		postclose `sim'

		preserve
		use `results',clear
		tempvar weight_mean
		 egen `weight_mean'=wtmean(spe),weight(N)
		qui sum `weight_mean',meanonly
		post `final_sim' (`dg') (r(mean))
		local mspe = string(round(r(mean),.0001))
		di as txt "polynomial degree=`dg' mspe=`mspe'"
		restore
	}
	
	
	postclose `final_sim'
}
ereturn clear							   
preserve
if("`method'"=="ife"){
	use `final_results',clear
	mkmat r mspe, matrix(CV)
	ereturn matrix CV = CV
	sort mspe
	qui sum mspe in 1 ,meanonly
	local min_mspe=r(mean)
	ereturn scalar min_mspe=`min_mspe'
	qui sum r in 1 ,meanonly
	local optimal_r=r(mean)
	ereturn scalar optimal_parameter=`optimal_r'
	if(`optimal_r'==0){
		ereturn scalar method=0
	}
	else{
		ereturn scalar method=1
	}
	di as res "optimal r=`optimal_r' in fe/ife model"
	restore	
}
if("`method'"=="mc"){
	use `final_results2',clear
	mkmat lambda mspe, matrix(CV)
	ereturn matrix CV = CV
	sort mspe
	qui sum mspe in 1,meanonly
	local min_mspe=r(mean)
	ereturn scalar min_mspe=`min_mspe'
	qui sum lambda in 1,meanonly
	local optimal_lambda=string(round(r(mean),.0001))
	ereturn scalar optimal_parameter=`optimal_lambda'

	qui sum lambda_norm in 1,meanonly
	local optimal_lambda_norm =string(round(r(mean),.0001))

	ereturn scalar method = 2
	di as res "optimal lambda=`optimal_lambda', lambda.norm=`optimal_lambda_norm' in mc model"
	restore
}
if("`method'"=="bspline"){
	use `final_results',clear
	mkmat nknots mspe, matrix(CV)
	ereturn matrix CV = CV
	sort mspe
	qui sum mspe in 1,meanonly
	local min_mspe=r(mean)
	ereturn scalar min_mspe=`min_mspe'
	qui sum nknots in 1,meanonly
	local optimal_nknots=string(round(r(mean),.0001))
	ereturn scalar optimal_parameter=`optimal_nknots'
	ereturn scalar method = 3
	di as res "optimal nknots=`optimal_nknots' in bspline model"
	restore
}
if("`method'"=="polynomial"){
	use `final_results',clear
	mkmat degree mspe, matrix(CV)
	ereturn matrix CV = CV
	sort mspe
	qui sum mspe in 1,meanonly
	local min_mspe=r(mean)
	ereturn scalar min_mspe=`min_mspe'
	qui sum degree in 1,meanonly
	local optimal_degree=r(mean)
	ereturn scalar optimal_parameter=`optimal_degree'
	ereturn scalar method = 4
	di as res "optimal degree=`optimal_degree' in polynomial model"
	restore
}
if("`method'"=="both"){
	use `final_results',clear
	mkmat r mspe, matrix(CV)
	matrix CV1 = CV
	sort mspe
	qui sum mspe in 1,meanonly
	local min_mspe1=r(mean)
	qui sum r in 1,meanonly
	local optimal_parameter1=r(mean)
	
	use `final_results2',clear
	mkmat lambda mspe, matrix(CV)
	matrix CV2 = CV
	sort mspe
	qui sum mspe in 1,meanonly
	local min_mspe2=r(mean)
	qui sum lambda in 1,meanonly
	local optimal_parameter2=r(mean)
	qui sum lambda_norm in 1,meanonly
	local optimal_parameter2_norm=r(mean)
	
	if(`min_mspe1'<=`min_mspe2'){
		di as res "Choose the fe/ife model, r = `optimal_parameter1'."
		ereturn matrix CV=CV1
		ereturn scalar min_mspe=`min_mspe1'
		if(`optimal_parameter1'==0){
			ereturn scalar method=0
		}
		else{
			ereturn scalar method=1
		}
		ereturn scalar optimal_parameter=`optimal_parameter1'
	}
	else{
		ereturn matrix CV=CV2
		ereturn scalar min_mspe=`min_mspe2'
		ereturn scalar optimal_parameter=`optimal_parameter2'
		local optimal_parameter2=string(round(`optimal_parameter2',.0001))
		local optimal_parameter2_norm=string(round(`optimal_parameter2_norm',.0001))
		di as res "Choose the mc model, lambda = `optimal_parameter2', lambda.norm = `optimal_parameter2_norm'."
		ereturn scalar method=2
	}
}
							   						   
end

program define placebo_Test,eclass
syntax varlist(min=1 max=1) , Treat(varlist min=1 max=1) ///
								  Unit(varlist min=1 max=1) ///
								  Time(varlist min=1 max=1) ///
								[ cov(varlist) ///
								  weight(varlist max=1) ///
								  force(string) ///
								  degree(integer 3) ///
								  nknots(integer 3) ///
								  r(integer 3) ///
								  lambda(real 0.5) ///
								  method(string) ///
								  alpha(real 0.05) ///
								  vartype(string) ///
								  nboots(integer 200) ///
								  tol(real 1e-4) ///
								  maxiterations(integer 5000) ///
								  seed(integer 12345678) ///
								  minT0(integer 5) ///
								  maxmissing(integer 0) ///
								  preperiod(integer -999999) ///
								  offperiod(integer 999999) ///
								  Proportion(real 0.3) ///
								  placeboperiod(integer 3) ///
								  Xlabel(string) ///
								  Ylabel(string) ///
								  Title(string) ///
								  Saving(string) ///
								]

set trace off
di as txt "{hline}"	
di as txt "Placebo Test..."
tempvar newid newtime touse
qui gen `touse'=1
qui  egen `newid'=group(`unit')
qui  egen `newtime'=group(`time')
sort `newid' `newtime'
local lag=`placeboperiod'

//generate placebotreat
tempvar ptreat treat_rev
qui mata:p_treat("`treat'","`newid'", "`newtime'","`ptreat'",`lag')
//generate true s
tempvar true_s true_s_off
qui gen `treat_rev' = 1 - `treat'
qui mata:gen_s("`treat'","`newid'", "`newtime'","`true_s'")
qui mata:gen_s("`treat_rev'","`newid'", "`newtime'","`true_s_off'")
//generate placebo s
tempvar placebo_t_on
qui mata:gen_s("`ptreat'","`newid'", "`newtime'","`placebo_t_on'")
//generate pllacebo s_off
tempvar ptreat_rev placebo_t_off
qui gen `ptreat_rev' = 1 - `ptreat'
qui mata:gen_s("`ptreat_rev'","`newid'", "`newtime'","`placebo_t_off'")

/* Threshold of touse */
tempvar ftreat notreatnum notmissingy missingy
gen `ftreat'=0
qui replace `ftreat'=1 if `ptreat'==0 & `varlist'!=. //the number of training data
bys `newid': egen `notreatnum'=sum(`ftreat')
qui replace `touse'=0 if `notreatnum'<`minT0'

if(`maxmissing'>0){
	bys `newid': egen `notmissingy'=count(`varlist')
	qui sum `newtime',meanonly
	local TT = r(max)-r(min)
	qui gen `missingy'=`TT'-`notmissingy'
	qui replace `touse'=0 if `missingy'>`maxmissing'
}
qui drop if `touse'==0
/* Threshold of touse END*/

preserve
tempvar counter TE S ATTs ATTs_off

qui fect_counter `counter' `TE' `ATTs' `ATTs_off',outcome(`varlist') treat(`ptreat') unit(`newid') time(`newtime') ///
t_on(`placebo_t_on') t_off(`placebo_t_off') t_on_original(`true_s') t_off_original(`true_s_off') ///
method("`method'") force("`force'") cov(`cov') nknots(`nknots') ///
degree(`degree') r(`r') tol(`tol')  ///
maxiterations(`maxiterations') lambda(`lambda')

qui sum `TE' [iweight = `weight'] if `true_s'<=0 & `true_s'>-`lag' & `touse'==1 ,meanonly
local PlaceboATT=r(mean)

forvalue s=`preperiod'/`offperiod' {

	qui sum `TE' [iweight = `weight'] if `true_s'==`s' & `touse'==1 ,meanonly
	
	if r(N)==0{
		di in red "No observations for s=`s'"
		local s_index=`s'-`preperiod'
		local att`s_index' = .
		local N`s_index'=0
	}
	else{
		local s_index=`s'-`preperiod'
		local att`s_index' = r(mean)
		local N`s_index'=r(N)
	}
	
}
restore
 
tempname sim sim2
tempfile results results2
qui postfile `sim' ATT using `results',replace
qui postfile `sim2' s ATTs using `results2',replace



if("`vartype'"=="bootstrap"){
	di as txt "Bootstrapping..."
	forvalue i=1/`nboots'{
		preserve
		tempvar counter TE S ATTs  ATTs_off
		bsample, cluster(`newid') idcluster(boot_id) 
		drop `newid'
		rename boot_id `newid'
			
		capture qui fect_counter `counter' `TE' `ATTs' `ATTs_off',outcome(`varlist') treat(`ptreat') unit(`newid') time(`newtime') ///
		t_on(`placebo_t_on') t_off(`placebo_t_off') t_on_original(`true_s') t_off_original(`true_s_off') ///
		method("`method'") force("`force'") cov(`cov') nknots(`nknots') degree(`degree') r(`r') tol(`tol') ///
		maxiterations(`maxiterations') lambda(`lambda')
			
		if _rc!=0 {
			restore
			continue
		}
			
		/* ATT */
		qui sum `TE' [iweight = `weight'] if `true_s'<=0 & `true_s'>-`lag' & `touse'==1,meanonly
		if r(mean)!=.{
			post `sim' (r(mean))
		}
			
		/* ATTs */
		forvalue si=`preperiod'/`offperiod' {
			qui sum `ATTs' [iweight = `weight'] if `true_s'==`si' & `touse'==1,meanonly
			if r(mean)!=.{
				post `sim2' (`si') (r(mean))
			}
		}
			
		if mod(`i',100)==0 {
			di as txt "ATT Estimation: Already Bootstrapped `i' Times"
		}
		restore
		}
		postclose `sim'
		postclose `sim2'
}
	
if("`vartype'"=="jackknife"){
	//di as txt "{hline}"	
	di as txt "Jackknifing..."
	tempvar order order_mean order_rank
	bys `newid': gen `order'=runiform()
	bys `newid':  egen `order_mean'=mean(`order')
	sort `order_mean'
	 egen `order_rank'=group(`order_mean')
	qui sum `newid',meanonly
	local max_id=r(max)
	forvalue i=1/`max_id'{
		preserve
		tempvar counter TE S ATTs ATTs_off
		qui drop if `order_rank'==`i'
			
		capture qui fect_counter `counter' `TE' `ATTs' `ATTs_off',outcome(`varlist') treat(`ptreat') unit(`newid') time(`newtime') ///
		t_on(`placebo_t_on') t_off(`placebo_t_off') t_on_original(`true_s') t_off_original(`true_s_off')  ///
		method("`method'") force("`force'") cov(`cov') nknots(`nknots') degree(`degree') r(`r') tol(`tol') ///
		maxiterations(`maxiterations') lambda(`lambda')
			
		if _rc!=0 {
			restore
			continue
		}
			
		/* ATT */
		qui sum `TE' [iweight = `weight'] if `true_s'<=0 & `true_s'>-`lag' & `touse'==1,meanonly
		if r(mean)!=.{
			post `sim' (r(mean))
		}
			
		/* ATTs */
		forvalue si=`preperiod'/`offperiod' {
			qui sum `ATTs' [iweight = `weight'] if `true_s'==`si' & `touse'==1,meanonly
			if r(mean)!=.{
				post `sim2' (`si') (r(mean))
			}
		}
			
		//di as txt "Jackknife: Round `i' Finished"
		restore
		}
		postclose `sim'
		postclose `sim2'
}
	
	
//Placebo_ATT
if("`vartype'"=="bootstrap"){
	preserve
	use `results',clear
	tempname ci placebo_ci_att
	qui sum ATT
	local Placebo_ATTsd=r(sd) //to save
	local twosidel=100*`alpha'/2
	local twosideu=100*(1-`alpha'/2)
	local twoside=100-`alpha'*100
	
	qui _pctile ATT,p(`twosidel' ,`twosideu')
	matrix `placebo_ci_att'=(r(r1),r(r2)) //to save
	matrix rownames `placebo_ci_att' = "Confidence Interval" 
	matrix colnames `placebo_ci_att' = "Lower Bound" "Upper Bound"
	if r(r1)>0 | r(r2)<0 {
		di "Placebo Test Fail"
	}
	else {
		di "Placebo Test Pass"
	}

	local placebo_z_score = `PlaceboATT'/`Placebo_ATTsd'
	local placebo_p_value = 2*(1-normal(abs(`placebo_z_score')))
	
	local placebo_pvalue=round(`placebo_p_value',0.001)
	local placebo_pvalueplot=string(`placebo_pvalue')
	

	//ATTs
	tempname sim_plot
	tempfile results_plot
	qui postfile `sim_plot' s atts att_lb att_ub s_N using `results_plot',replace
	
	qui use `results2',clear

	forvalue s=`preperiod'/`offperiod' {
		qui sum ATTs if s==`s'
		local attsd_temp=r(sd)
		qui _pctile ATTs if s==`s',p(`twosidel' ,`twosideu')
		matrix `ci'=(r(r1),r(r2))
		local s_index=`s'-`preperiod'
		//post `sim_plot' (`s') (`att`s_index'') (`ci'[1,1]) (`ci'[1,2]) (`N`s_index'')
		
		if("`vartype'"=="bootstrap"){
			post `sim_plot' (`s') (`att`s_index'') (`ci'[1,1]) (`ci'[1,2]) (`N`s_index'')
		}
		if("`vartype'"=="jackknife"){
			local tvalue = invttail(`max_id'-1,`alpha'/2)
			local ub_jack = `att`s_index''+`tvalue'*`attsd_temp'
			local lb_jack = `att`s_index''-`tvalue'*`attsd_temp'
			post `sim_plot' (`s') (`att`s_index'') (`lb_jack') (`ub_jack') (`N`s_index'')
		}
	}
	postclose `sim_plot'
	restore
}


if("`vartype'"=="jackknife"){
	preserve
	use `results',clear
	tempname ci placebo_ci_att
	qui sum ATT
	local Placebo_ATTsd=r(sd)*sqrt(`max_id') //to save
	
	local tvalue = invttail(`max_id'-1,`alpha'/2)
	local ub_jack_p = `PlaceboATT'+`tvalue'*`Placebo_ATTsd'
	local lb_jack_p = `PlaceboATT'-`tvalue'*`Placebo_ATTsd'
		
	matrix `placebo_ci_att'=(`lb_jack_p',`ub_jack_p') //to save
	matrix rownames `placebo_ci_att' = "Confidence Interval" 
	matrix colnames `placebo_ci_att' = "Lower Bound" "Upper Bound"
	

	//p-value
	local placebo_z_score = `PlaceboATT'/`Placebo_ATTsd'
	local placebo_p_value = 2*(1-normal(abs(`placebo_z_score')))
	
	local placebo_pvalue=round(`placebo_p_value',0.001)
	local placebo_pvalueplot=string(`placebo_pvalue')

	//ATTs
	tempname sim_plot
	tempfile results_plot
	qui postfile `sim_plot' s atts att_lb att_ub s_N using `results_plot',replace
	
	qui use `results2',clear

	forvalue s=`preperiod'/`offperiod' {
		qui sum ATTs if s==`s'
		local attsd_temp=r(sd)*sqrt(`max_id')
		//qui _pctile ATTs if s==`s',p(`twosidel' ,`twosideu')
		//matrix `ci'=(r(r1),r(r2))
		local s_index=`s'-`preperiod'
		//post `sim_plot' (`s') (`att`s_index'') (`ci'[1,1]) (`ci'[1,2]) (`N`s_index'')
		
		if("`vartype'"=="jackknife"){
			local tvalue = invttail(`max_id'-1,`alpha'/2)
			local ub_jack = `att`s_index''+`tvalue'*`attsd_temp'
			local lb_jack = `att`s_index''-`tvalue'*`attsd_temp'
			post `sim_plot' (`s') (`att`s_index'') (`lb_jack') (`ub_jack') (`N`s_index'')
		}
	}
	postclose `sim_plot'
	restore
}



// Draw ATTs plot
preserve
qui use `results_plot',clear

qui  egen max_N=max(s_N)
qui sum max_N,meanonly
local max_N=r(mean)
/*NEW update with regard to proportion*/
local N_cutoff = `max_N'*`proportion'

qui drop if s_N<`N_cutoff'
qui drop if s<`preperiod'
qui drop if s>`offperiod'

//update preperiod and offperiod
qui sum s, meanonly
local preperiod = r(min)
local offperiod = r(max)
/*NEW End*/

qui tsset s
local CIlevel=100*(1-`alpha')
qui  egen max_atts=max(att_ub)

local axis2_up=4*`max_N'
qui sum max_atts,meanonly
local max_atts=1.2*r(mean)
if `max_atts'<0 {
	local max_atts=0
}
qui  egen min_atts=min(att_lb)
qui sum min_atts,meanonly
local min_atts=1.2*r(mean)
if `min_atts'>0 {
	local min_atts=0
}
qui gen y0=0
local placebo_line=`lag'-1+0.5

if(abs(`max_atts')>abs(`min_atts')){
		local min_atts=-`max_atts'*0.8
	}
	else{
		local max_atts=-`min_atts'*0.8
}

/* Title */
local xlabel = cond("`xlabel'"=="","Time relative to the Treatment","`xlabel'")
local ylabel = cond("`ylabel'"=="","Average Treatment Effect","`ylabel'")
local title = cond("`title'"=="","Placebo Test","`title'")
local neg_lag=-(`lag'-1)


/* NEW Changes LSJ */
twoway (rcap att_lb att_ub s, color(black) lcolor(black) yaxis(1)) ///
(scatter atts s, mcolor(black) yaxis(1) msize(2pt)) ///
(rcap att_lb att_ub s if s>-`lag' & s<=0, color(blue) lcolor(blue) yaxis(1)) ///
(scatter atts s if s>-`lag' & s<=0,mcolor(blue) msymbol(O) yaxis(1) msize(small)) ///
	(bar s_N s,yaxis(2) barwidth(0.5) color(gray%50)), ///
	 xline(0.5 -`placebo_line', lcolor(gs6%50) lpattern(dash)) ///
	 yline(0, lcolor(gs6%50) lpattern(dash)) ///
	 yscale(axis(2) r(0 `axis2_up')) ///
	 yscale(axis(1) r(`min_atts' `max_atts')) ///
	 xscale(noextend) ///
	 xlabel(`preperiod' 0 (10) `offperiod', labsize(*0.75)) ///
	 xtick(`preperiod'(1)`offperiod') ///
	 ylabel(0 `max_N',axis(2) ) ///
	 ytitle(`ylabel',axis(1)) ///
	 ytitle("Num of observations",axis(2)) ///
	 xtitle(`xlabel') ///
	 title(`title') ///
	 text(`max_atts' `offperiod' "Placebo Test p-value: `placebo_pvalueplot'",place(sw)) ///
	 scheme(s2mono) ///
	 graphr(fcolor(white)) ///
	 plotregion(lcolor(black) lwidth(thin) margin(zero) ) ///
	 note("`CIlevel'% Confidence Interval") ///
legend(order( 1 "ATT(`CIlevel'% CI)" 2 "Placebo Region" 3 "ATT")) 


if("`saving'"!=""){
   cap graph export `saving',replace
} 
restore


ereturn scalar placebo_pvalue=`placebo_pvalue'
ereturn scalar placebo_ATT_mean=`PlaceboATT'
ereturn scalar placebo_ATT_sd=`Placebo_ATTsd'
ereturn matrix placebo_CI=`placebo_ci_att'


end

/*NEW LSJ Carryover Effect Test*/
program define carryover_Test,eclass
syntax varlist(min=1 max=1) , Treat(varlist min=1 max=1) ///
								  Unit(varlist min=1 max=1) ///
								  Time(varlist min=1 max=1) ///
								[ cov(varlist) ///
								  weight(varlist max=1) ///
								  force(string) ///
								  degree(integer 3) ///
								  nknots(integer 3) ///
								  r(integer 3) ///
								  lambda(real 0.5) ///
								  method(string) ///
								  alpha(real 0.05) ///
								  vartype(string) ///
								  nboots(integer 200) ///
								  tol(real 1e-4) ///
								  maxiterations(integer 5000) ///
								  seed(integer 12345678) ///
								  minT0(integer 5) ///
								  maxmissing(integer 0) ///
								  preperiod(integer -999999) ///
								  offperiod(integer 999999) ///
								  Proportion(real 0.3) ///
								  carryoverperiod(integer 3) ///
								  Xlabel(string) ///
								  Ylabel(string) ///
								  Title(string) ///
								  Saving(string) ///
								]

set trace off
di as txt "{hline}"	
di as txt "Carryover Test..."
tempvar newid newtime touse
qui gen `touse'=1
qui  egen `newid'=group(`unit')
qui  egen `newtime'=group(`time')
sort `newid' `newtime'
local lead=`carryoverperiod'

//generate carryovertreat
tempvar ctreat
qui mata:c_treat("`treat'","`newid'", "`newtime'","`ctreat'",`lead')
//generate true s
tempvar true_s
qui mata:gen_s("`treat'","`newid'", "`newtime'","`true_s'")
//generate true reverse s
tempvar treat_rev true_off_s
qui gen `treat_rev' = 1 - `treat'
qui mata:gen_s("`treat_rev'","`newid'", "`newtime'","`true_off_s'")
//generate carryover s
tempvar carryover_t_on
qui mata:gen_s("`ctreat'","`newid'", "`newtime'","`carryover_t_on'")
//generate carryover s_off
tempvar ctreat_rev carryover_t_off
qui gen `ctreat_rev' = 1 - `ctreat'
qui mata:gen_s("`ctreat_rev'","`newid'", "`newtime'","`carryover_t_off'")

/* Threshold of touse */
tempvar ftreat notreatnum notmissingy missingy
gen `ftreat'=0
qui replace `ftreat'=1 if `ctreat'==0 & `varlist'!=. //the number of training data
bys `newid': egen `notreatnum'=sum(`ftreat')
qui replace `touse'=0 if `notreatnum'<`minT0'

if(`maxmissing'>0){
	bys `newid': egen `notmissingy'=count(`varlist')
	qui sum `newtime',meanonly
	local TT = r(max)-r(min)
	qui gen `missingy'=`TT'-`notmissingy'
	qui replace `touse'=0 if `missingy'>`maxmissing'
}
qui drop if `touse'==0
/* Threshold of touse END*/

preserve
tempvar counter TE S ATTs ATTs_off

qui fect_counter `counter' `TE' `ATTs' `ATTs_off',outcome(`varlist') treat(`ctreat') unit(`newid') time(`newtime') ///
t_on(`carryover_t_on') t_off(`carryover_t_off') t_on_original(`true_s') t_off_original(`true_off_s') ///
method("`method'") force("`force'") cov(`cov') nknots(`nknots') ///
degree(`degree') r(`r') tol(`tol')  ///
maxiterations(`maxiterations') lambda(`lambda')
/*QUESTION: lead 3 periods, including 0?*/
qui sum `TE' [iweight = `weight'] if `true_off_s'>0 & `true_off_s'<=`lead' & `touse'==1 ,meanonly
local CarryoverATT=r(mean)

forvalue s=`preperiod'/`offperiod' {

	qui sum `TE' [iweight = `weight'] if `true_off_s'==`s' & `touse'==1 ,meanonly
	
	if r(N)==0{
		di in red "No observations for s=`s'"
		local s_index=`s'-`preperiod'
		local att_off`s_index' = .
		local N_off`s_index'=0
	}
	else{
		local s_index=`s'-`preperiod'
		local att_off`s_index' = r(mean)
		local N_off`s_index'=r(N)
	}
	
}
restore
 
tempname sim sim2_off
tempfile results results2_off
qui postfile `sim' ATT using `results',replace
qui postfile `sim2_off' s ATTs using `results2_off',replace
// sim-results: ATT only in 0-lead periods
// sim2_off-results2_off: ATT for all the periods

if("`vartype'"=="bootstrap"){
	di as txt "Bootstrapping..."
	forvalue i=1/`nboots'{
		preserve
		tempvar counter TE S ATTs  ATTs_off
		bsample, cluster(`newid') idcluster(boot_id) 
		drop `newid'
		rename boot_id `newid'
			
		capture qui fect_counter `counter' `TE' `ATTs' `ATTs_off',outcome(`varlist') treat(`ctreat') unit(`newid') time(`newtime') ///
		t_on(`carryover_t_on') t_off(`carryover_t_off') t_on_original(`true_s') t_off_original(`true_off_s')  ///
		method("`method'") force("`force'") cov(`cov') nknots(`nknots') degree(`degree') r(`r') tol(`tol') ///
		maxiterations(`maxiterations') lambda(`lambda')
			
		if _rc!=0 {
			restore
			continue
		}
			
		/* ATT */
		qui sum `TE' [iweight = `weight'] if `true_off_s'>0 & `true_off_s'<=`lead' & `touse'==1,meanonly
		if r(mean)!=.{
			post `sim' (r(mean))
		}
			
		forvalue si=`preperiod'/`offperiod' {
			qui sum `ATTs_off' [iweight = `weight'] if `true_off_s'==`si' & `touse'==1,meanonly
			if r(mean)!=.{
				post `sim2_off' (`si') (r(mean))
			}
		}
			
		if mod(`i',100)==0 {
			di as txt "ATT Estimation: Already Bootstrapped `i' Times"
		}
		restore
		}
		postclose `sim'
		postclose `sim2_off'
}
	
if("`vartype'"=="jackknife"){
	//di as txt "{hline}"	
	di as txt "Jackknifing..."
	tempvar order order_mean order_rank
	bys `newid': gen `order'=runiform()
	bys `newid':  egen `order_mean'=mean(`order')
	sort `order_mean'
	 egen `order_rank'=group(`order_mean')
	qui sum `newid',meanonly
	local max_id=r(max)
	forvalue i=1/`max_id'{
		preserve
		tempvar counter TE S ATTs ATTs_off
		qui drop if `order_rank'==`i'
			
		capture qui fect_counter `counter' `TE' `ATTs' `ATTs_off',outcome(`varlist') treat(`ctreat') unit(`newid') time(`newtime') ///
		t_on(`carryover_t_on') t_off(`carryover_t_off') t_on_original(`true_s') t_off_original(`true_off_s')   ///
		method("`method'") force("`force'") cov(`cov') nknots(`nknots') degree(`degree') r(`r') tol(`tol') ///
		maxiterations(`maxiterations') lambda(`lambda')
			
		if _rc!=0 {
			restore
			continue
		}
			
		/* ATT */
		qui sum `TE' [iweight = `weight'] if `true_off_s'>0 & `true_off_s'<=`lead' & `touse'==1,meanonly
		if r(mean)!=.{
			post `sim' (r(mean))
		}
			
		/* ATTs */
		forvalue si=`preperiod'/`offperiod' {
			qui sum `ATTs_off' [iweight = `weight'] if `true_off_s'==`si' & `touse'==1,meanonly
			if r(mean)!=.{
				post `sim2_off' (`si') (r(mean))
			}
		}
			
		//di as txt "Jackknife: Round `i' Finished"
		restore
		}
		postclose `sim'
		postclose `sim2_off'
}
	
	
//Carryover_ATT
if("`vartype'"=="bootstrap"){
	preserve
	use `results',clear
	tempname ci carryover_ci_att
	qui sum ATT
	local Carryover_ATTsd=r(sd) //to save
	local twosidel=100*`alpha'/2 // twoside test
	local twosideu=100*(1-`alpha'/2)
	local twoside=100-`alpha'*100
	
	qui _pctile ATT,p(`twosidel' ,`twosideu')
	matrix `carryover_ci_att'=(r(r1),r(r2)) //to save
	matrix rownames `carryover_ci_att' = "Confidence Interval" 
	matrix colnames `carryover_ci_att' = "Lower Bound" "Upper Bound"
	if r(r1)>0 | r(r2)<0 {
		di "Carryover Test Fail"
	}
	else {
		di "Carryover Test Pass"
	}

	local carryover_z_score = `CarryoverATT'/`Carryover_ATTsd'
	local carryover_p_value = 2*(1-normal(abs(`carryover_z_score')))
	
	local carryover_pvalue=round(`carryover_p_value',0.001)
	local carryover_pvalueplot=string(`carryover_pvalue')
	

	//ATTs
	tempname sim_plot_off
	tempfile results_plot_off 
	//next line no attsd?
	qui postfile `sim_plot_off' s atts att_lb att_ub s_N using `results_plot_off',replace
	
	qui use `results2_off',clear

	forvalue s=`preperiod'/`offperiod' {
		qui sum ATTs if s==`s'
		local attsd_temp=r(sd)
		qui _pctile ATTs if s==`s',p(`twosidel' ,`twosideu')
		matrix `ci'=(r(r1),r(r2))
		local s_index=`s'-`preperiod'
		//post `sim_plot' (`s') (`att`s_index'') (`ci'[1,1]) (`ci'[1,2]) (`N`s_index'')
		// must be something wrong here
		if("`vartype'"=="bootstrap"){
			post `sim_plot_off' (`s') (`att_off`s_index'') (`ci'[1,1]) (`ci'[1,2]) (`N_off`s_index'')
		}
	}
	postclose `sim_plot_off'
	restore
}


if("`vartype'"=="jackknife"){
	preserve
	use `results',clear
	tempname ci carryover_ci_att
	qui sum ATT
	local Carryover_ATTsd=r(sd)*sqrt(`max_id') //to save
	
	local tvalue = invttail(`max_id'-1,`alpha'/2)
	local ub_jack_c = `CarryoverATT'+`tvalue'*`Carryover_ATTsd'
	local lb_jack_c = `CarryoverATT'-`tvalue'*`Carryover_ATTsd'
		
	matrix `carryover_ci_att'=(`lb_jack_c',`ub_jack_c') //to save
	matrix rownames `carryover_ci_att' = "Confidence Interval" 
	matrix colnames `carryover_ci_att' = "Lower Bound" "Upper Bound"
	

	//p-value
	local carryover_z_score = `CarryoverATT'/`Carryover_ATTsd'
	local carryover_p_value = 2*(1-normal(abs(`carryover_z_score')))
	
	local carryover_pvalue=round(`carryover_p_value',0.001)
	local carryover_pvalueplot=string(`carryover_pvalue')

	//ATTs
	tempname sim_plot_off
	tempfile results_plot_off
	qui postfile `sim_plot_off' s atts att_lb att_ub s_N using `results_plot_off',replace
	
	qui use `results2_off',clear

	forvalue s=`preperiod'/`offperiod' {
		qui sum ATTs if s==`s'
		local attsd_temp=r(sd)*sqrt(`max_id')
		//qui _pctile ATTs if s==`s',p(`twosidel' ,`twosideu')
		//matrix `ci'=(r(r1),r(r2))
		local s_index=`s'-`preperiod'
		//post `sim_plot' (`s') (`att`s_index'') (`ci'[1,1]) (`ci'[1,2]) (`N`s_index'')
		
		if("`vartype'"=="jackknife"){
			local tvalue = invttail(`max_id'-1,`alpha'/2)
			local ub_jack = `att_off`s_index''+`tvalue'*`attsd_temp'
			local lb_jack = `att_off`s_index''-`tvalue'*`attsd_temp'
			post `sim_plot_off' (`s') (`att_off`s_index'') (`lb_jack') (`ub_jack') (`N_off`s_index'')
		}
	}
	postclose `sim_plot_off'
	restore
}

// Draw ATTs plot
preserve
qui use `results_plot_off',clear

qui  egen max_N=max(s_N)
qui sum max_N,meanonly
local max_N=r(mean)
/*NEW update with regard to proportion*/
local N_cutoff = `max_N'*`proportion'

qui drop if s_N<`N_cutoff'
qui drop if s<`preperiod'
qui drop if s>`offperiod'

//update preperiod and offperiod
qui sum s, meanonly
local preperiod = r(min)
local offperiod = r(max)
/*NEW End*/

qui tsset s
local CIlevel=100*(1-`alpha')
qui  egen max_atts=max(att_ub)

local axis2_up=4*`max_N'
qui sum max_atts,meanonly
local max_atts=1.2*r(mean)
if `max_atts'<0 {
	local max_atts=0
}
qui  egen min_atts=min(att_lb)
qui sum min_atts,meanonly
local min_atts=1.2*r(mean)
if `min_atts'>0 {
	local min_atts=0
}
qui gen y0=0
local carryover_line=`lead'+0.5

if(abs(`max_atts')>abs(`min_atts')){
		local min_atts=-`max_atts'*0.8
	}
	else{
		local max_atts=-`min_atts'*0.8
}

/* Title */
local xlabel = cond("`xlabel'"=="","Time relative to the Exit of Treatment","`xlabel'")
local ylabel = cond("`ylabel'"=="","Average Treatment Effect","`ylabel'")
local title = cond("`title'"=="","Carryover Test","`title'")
local neg_lead=`lead'+1
/*MEW LSJ*/
twoway (rcap att_lb att_ub s, color(black) lcolor(black) yaxis(1)) ///
(rcap att_lb att_ub s if s<=`lead' & s>0, color(red) lcolor(red) yaxis(1)) ///
(scatter atts s, mcolor(black) yaxis(1) msize(2pt)) ///
(scatter atts s if s<=`lead' & s>0,mcolor(red) msymbol(O) yaxis(1) msize(small)) ///
(bar s_N s,yaxis(2) barwidth(0.5) color(gray%50)), ///
xline(0.5 `carryover_line', lcolor(gs6%50) lpattern(dash)) ///
yline(0, lcolor(gs6%50) lpattern(dash)) ///
ysc(axis(1) r(`min_atts' `max_atts')) ///
ysc(axis(2) r(0 `axis2_up')) ///
ylabel(0 `max_N',axis(2)) ///
xlabel(`preperiod' 0 (10) `neg_lead' `offperiod', labsize(*0.75)) ///
xtick(`preperiod'(1)`offperiod') ///
ytitle(`ylabel',axis(1)) ///
ytitle("Num of observations",axis(2)) ///
xtitle(`xlabel') ///
title(`title') ///
text(`max_atts' `offperiod' "Carryover Test p-value: `carryover_pvalueplot'",place(sw)) ///
scheme(s2mono) ///
graphr(fcolor(white)) ///
plotregion(lcolor(black) lwidth(thin) margin(zero) ) ///
note("`CIlevel'% Confidence Interval") ///
legend(order( 1 "ATT(`CIlevel'% CI)" 2 "Carryover Region" 3 "ATT")) 


if("`saving'"!=""){
   cap graph export `saving',replace
} 
restore


ereturn scalar carryover_pvalue=`carryover_pvalue'
ereturn scalar carryover_ATT_mean=`CarryoverATT'
ereturn scalar carryover_ATT_sd=`Carryover_ATTsd'
ereturn matrix carryover_CI=`carryover_ci_att'


end
/*LSJ END*/


program define permutation,eclass
syntax varlist(min=1 max=1) , Treat(varlist min=1 max=1) ///
								  Unit(varlist min=1 max=1) ///
								  Time(varlist min=1 max=1) ///
								  t_on(varlist min=1 max=1) ///
								  t_off(varlist min=1 max=1) ///
								[ cov(varlist) ///
								  weight(varlist max=1) ///
								  force(string) ///
								  degree(integer 3) ///
								  nknots(integer 3) ///
								  r(integer 3) ///
								  lambda(real 0.5) ///
								  method(string) ///
								  nboots(integer 200) ///
								  tol(real 1e-4) ///
								  maxiterations(integer 5000) ///
								]
								
set trace off
di as txt "{hline}"	
di as txt "Permutation Test..."
tempvar newid newtime
qui  egen `newid'=group(`unit')
qui  egen `newtime'=group(`time')
sort `newid' `newtime'

/* estimate ATT */
preserve
tempvar counter TE S ATTs ATTs_off

qui fect_counter `counter' `TE' `ATTs' `ATTs_off',outcome(`varlist') treat(`treat') unit(`newid') time(`newtime') ///
t_on(`t_on') t_off(`t_off') ///
method("`method'") force("`force'") cov(`cov') nknots(`nknots') degree(`degree') r(`r') tol(`tol') ///
maxiterations(`maxiterations') lambda(`lambda')

sort `newid' `newtime'
qui sum `TE' [iweight = `weight'] if `treat'==1,meanonly
local ATT=r(mean) //to save

restore

/* Simulation */
tempname sim_permu
tempfile results_permu
qui postfile `sim_permu' ATT using `results_permu',replace

forvalue i=1/`nboots'{
	preserve
	tempvar permutreat
	qui mata:permute("`treat'","`newid'","`newtime'","`permutreat'")
	
	tempvar counter TE S ATTs ATTs_off
	qui fect_counter `counter' `TE' `ATTs' `ATTs_off',outcome(`varlist') treat(`permutreat') unit(`newid') time(`newtime') ///
	t_on(`t_on') t_off(`t_off')  ///
	method("`method'") force("`force'") cov(`cov') nknots(`nknots') degree(`degree') r(`r') tol(`tol') ///
	maxiterations(`maxiterations') lambda(`lambda')
	
	qui sum `TE' [iweight = `weight'] if `permutreat'==1,meanonly
	if r(mean)!=.{
		post `sim_permu' (r(mean))
	}
	
	if mod(`i',100)==0 {
			di as txt "Permutation Test: Already Simulated `i' Times"
	}
	restore
}
postclose `sim_permu'

preserve
qui use `results_permu',clear
qui drop if ATT<`ATT'
qui sum ATT,meanonly
local pvalue=r(N)/`nboots'
ereturn scalar permutation_pvalue=`pvalue'
local pvalue_print=string(round(`pvalue',.001))
di as res "Permutation Test: p value=`pvalue_print'"
restore
end

/* program define _gwtmean
	version 3.0
	local varlist "req new max(1)"
	local exp "req nopre"
	local if "opt"
	local in "opt"
	local options "by(string) weight(string)"
	parse "`*'"
	tempvar touse 
	if "`weight'" ~= "" {
		local weight "* (`weight')"
	}
	quietly {
		gen byte `touse'=1 `if' `in'
		sort `touse' `by'
		by `touse' `by': replace `varlist' = /*
			*/ sum((`exp')`weight')/sum(((`exp')!=.)`weight') if `touse'==1
		by `touse' `by': replace `varlist' = `varlist'[_N]
	}
end */

*********************************************************mata
mata:

void treat_fill(string missing_treat,string unit, string time)
{
	real matrix treat,panel_treat,x
	real scalar minid,maxid,mintime,maxtime,rownum,colnum,i,j
	st_view(treat=.,.,missing_treat)
	st_view(x=.,.,(unit,time))

	minid=min(x[,1])
	maxid=max(x[,1])
	mintime=min(x[,2])
	maxtime=max(x[,2])
	rownum=maxid-minid+1
	colnum=maxtime-mintime+1

	panel_treat=treat
	panel_treat=rowshape(panel_treat,rownum)

	for(i=1;i<=rownum;i++) {
		if(hasmissing(panel_treat[i,1])==1) {
			panel_treat[i,1]=0
		}
	
		for(j=2;j<=colnum;j++) {
			if(hasmissing(panel_treat[i,j])==1) {
				panel_treat[i,j]=panel_treat[i,j-1]
			}
		}
	}
	treat=colshape(panel_treat,1)
	st_store(.,missing_treat,treat)
}


void gen_s(string missing_treat,string unit, string time, string outputs)
{
	real matrix treat,panel_treat,x,panel_s
	real scalar minid,maxid,mintime,maxtime,rownum,colnum,i,j,state,start,cut,k
	st_view(treat=.,.,missing_treat)
	st_view(x=.,.,(unit,time))

	minid=min(x[,1])
	maxid=max(x[,1])
	mintime=min(x[,2])
	maxtime=max(x[,2])
	rownum=maxid-minid+1
	colnum=maxtime-mintime+1

	panel_treat=treat
	panel_treat=rowshape(panel_treat,rownum)
	panel_s=panel_treat

	for(i=1;i<=rownum;i++) {
		if(panel_treat[i,1]==0){
			state=0
			start=1
			cut=1
			panel_s[i,1]=.
		}
		else{
			state=1
			panel_s[i,1]=.
		}
		for(j=2;j<=colnum;j++) {
			if(panel_treat[i,j]==0 & state==1) {
				state=0
				start=j
				cut=j
				panel_s[i,j]=.
			}
			else if(panel_treat[i,j]==0 & state==0) {
				cut=cut+1
				panel_s[i,j]=.
			}
			else if(panel_treat[i,j]==1 & state==1) {
				panel_s[i,j]=1+panel_s[i,j-1]
			}
			else if(panel_treat[i,j]==1 & state==0){
				for(k=start;k<=cut;k++) {
					panel_s[i,k]=k-cut
				}
				panel_s[i,j]=1
				state=1
			}
		}
	}
	panel_s=colshape(panel_s,1)
	st_addvar("int", outputs)
	st_store(., outputs,panel_s)
}

/* NEW LSJ: Invert treatment: 1 to 0, 0 to 1 */
void inv_treat(string Treat,
			 string unit, 
			 string time, 
			 string outputs)
{
	real matrix treat,panel_treat,x,inverttreat
	real scalar minid,maxid,mintime,maxtime,rownum,colnum,i,j,k
	st_view(treat=.,.,Treat)
	st_view(x=.,.,(unit,time))

	minid=min(x[,1])
	maxid=max(x[,1])
	mintime=min(x[,2])
	maxtime=max(x[,2])
	rownum=maxid-minid+1
	colnum=(maxtime-mintime)+1

	panel_treat=treat
	panel_treat=rowshape(panel_treat,rownum)
    inverttreat=panel_treat

	for(i=1;i<=rownum;i++){
        for(j=colnum;j>=1;j--){
            if(panel_treat[i,j]==1){
                inverttreat[i,j]=0
            }
            else if(panel_treat[i,j]==0){
                inverttreat[i,j]=1
            }
        }
    }
	inverttreat=colshape(inverttreat,1) 
	st_addvar("int", outputs)
	st_store(., outputs,inverttreat)
}

void ife(string outcome, 
         string cov, 
         string treat,  
         string unit, 
         string time,
         string force, 
         real r, 
         real Tol, 
         real iter, 
         string newvarname, 
         string alpha, 
         string xi, 
         string mu)
{   
    real matrix D_long,D_panel
    real matrix I_long,I_panel
    real matrix X_long,Y_long,Y_panel,Y_predict,Y_predict_old,Y_predict_delta
    real scalar num_obs,p
    real matrix unit_long, time_long
    real scalar num_unit,num_time
    real matrix X_long_I, XX_inv, XY, beta_hat
    real matrix alpha_hat, xi_hat, mu_hat, ife_hat
    real matrix Y_hat_panel
    real matrix Y_hat_long_I, Y_hat_panel_I
    real matrix W
    real matrix W_all_mean, W_unit_mean, W_time_mean, W_hat
    real matrix I_NT
    real scalar trans,m,crit,i
    real matrix U, Vt, svd_val,target_svd,ife_fit

    st_view(D_long=.,.,treat)
	D_long = editmissing(D_long,0)
    I_long=(D_long:==0)
    st_view(X_long=.,.,(cov))
	X_long = editmissing(X_long,0)

    st_view(Y_long=.,.,outcome)
	Y_long = editmissing(Y_long,0)
	num_obs=rows(Y_long)
	p=cols(X_long)

    st_view(unit_long=.,.,unit)
	st_view(time_long=.,.,time)
	unit_long = editmissing(unit_long,0)
	time_long = editmissing(time_long,0)
	num_unit=rows(uniqrows(unit_long))
	num_time=rows(uniqrows(time_long))

    I_NT = J(num_unit,num_time,1)

    //get XX_inv (p×p dimension)
	if(p>0){
		X_long_I = X_long:*I_long //NT*p
		XX_inv = transposeonly(X_long_I)*X_long_I
		_invsym(XX_inv) //p*p
	}
	
    //reshape Y_long; unit_long; time_long
    Y_panel=rowshape(Y_long,num_unit) //N*T
    D_panel=rowshape(D_long,num_unit) //N*T
    I_panel=rowshape(I_long,num_unit) //N*T
	
    //initialization
	st_view(alpha_hat=.,.,alpha)
	alpha_hat = editmissing(alpha_hat,0)
    alpha_hat = rowshape(alpha_hat,num_unit) //N*T
	st_view(xi_hat=.,.,xi)
	xi_hat = editmissing(xi_hat,0)
    xi_hat = rowshape(xi_hat,num_unit) //N*T
	st_view(mu_hat=.,.,mu)	
	mu_hat = editmissing(mu_hat,0)
    mu_hat = rowshape(mu_hat,num_unit)  //N*T
	ife_hat = J(num_unit,num_time,0)  //N*T

    // EM Algorithm
    for(m=1;m<=iter;m++) {
        Y_hat_panel = Y_panel-mu_hat-alpha_hat-xi_hat-ife_hat
        Y_hat_panel_I = Y_hat_panel:*I_panel
        Y_hat_long_I = rowshape(Y_hat_panel_I,num_obs)
        if(p>0){
			XY = transposeonly(X_long_I)*Y_hat_long_I
			beta_hat = XX_inv*XY
		}
		
        //calculate W, a N*T matrix
        W=J(num_unit,num_time,0) //N*T
        if(p>0){
            W = (Y_panel - rowshape(X_long*beta_hat,num_unit)):*I_panel
        }else{
            W = Y_panel:*I_panel
        }
        W = W + (mu_hat  + alpha_hat + xi_hat + ife_hat):*D_panel
        
        //calculate W_mean; W_unit_mean; W_time_mean
        W_all_mean = J(num_unit,num_time,sum(W)/num_obs)
        W_unit_mean = J(num_unit,num_time,0)
        W_time_mean = J(num_unit,num_time,0)
        if(force=="two-way"|force=="unit"){ // mean across T
			W_unit_mean = rowsum(W):*I_NT/num_time
		}
		
        if(force=="two-way"|force=="time"){ // mean across N
			W_time_mean = I_NT:*colsum(W)/num_unit
		}
        if(force=="two-way"){
			W_hat = W - (W_unit_mean+W_time_mean) + W_all_mean
		}
		if(force=="unit"){
			W_hat = W - W_unit_mean
		}
		if(force=="time"){
			W_hat = W - W_time_mean
		}
		if(force=="none"){
			W_hat = W - W_all_mean
		}
		
        //svd
		if(rows(W_hat)<cols(W_hat)){
			trans=1
			W_hat = transposeonly(W_hat)
		}else{
			trans=0
		}
		W_hat = W_hat/num_obs
		svd(W_hat,U=.,svd_val=.,Vt=.)
		target_svd=J(length(svd_val),1,0)
		for(i=1;i<=r;i++){
			target_svd[i,]=1
		}
		svd_val=svd_val:*target_svd
		ife_fit = U[|1,1\.,r|]*diag(svd_val)[|1,1\r,r|]*Vt[|1,1\r,.|]
		if(trans==1){
			ife_fit=transposeonly(ife_fit)
		}
		
        //update
		ife_hat=ife_fit*num_obs
		if(force=="two-way"|force=="unit"){
			alpha_hat = W_unit_mean-W_all_mean
		}
		if(force=="two-way"|force=="time"){
			xi_hat = W_time_mean - W_all_mean
		}
		mu_hat = W_all_mean

        //predict
		Y_predict = mu_hat+alpha_hat+xi_hat+ife_hat
		if(p>0){
			Y_predict = Y_predict + rowshape(X_long*beta_hat,num_unit)
		}

        if(m==1){
			Y_predict_old = Y_predict
		}else{ 
			Y_predict_delta = (Y_predict-Y_predict_old):*I_panel
			crit = sqrt(sum(Y_predict_delta:*Y_predict_delta))/sqrt(sum(Y_predict:*Y_predict))
			if (crit<Tol){
				st_numscalar("r(stop)", 1)
				break
			} 
			st_numscalar("r(stop)", 0)
			Y_predict_old = Y_predict
		}
    }
	Y_predict = rowshape(Y_predict,num_obs)
    st_addvar("double", newvarname)
	st_store(., newvarname,Y_predict)

	if(p>0){
		st_matrix("r(coef)",transposeonly(beta_hat))
	}
	st_numscalar("r(cons)", sum(W_all_mean)/num_obs)
}


void mc(string outcome, 
        string cov, 
        string treat,  
        string unit, 
        string time,
        string force ,
        real lambda, 
        real Tol, 
        real iter, 
        string newvarname, 
        string alpha, 
        string xi, 
        string mu)
{    
    real matrix D_long,D_panel
    real matrix I_long,I_panel
    real matrix X_long,Y_long,Y_panel,Y_predict,Y_predict_old,Y_predict_delta
    real scalar num_obs,p
    real matrix unit_long, time_long
    real scalar num_unit,num_time
    real matrix X_long_I, XX_inv, XY, beta_hat
    real matrix alpha_hat, xi_hat, mu_hat, mc_hat
    real matrix Y_hat_panel
    real matrix Y_hat_long_I, Y_hat_panel_I
    real matrix W
    real matrix W_all_mean, W_unit_mean, W_time_mean, W_hat
    real matrix I_NT
    real scalar trans,m,crit
    real matrix U, Vt, svd_val,mc_fit

    st_view(D_long=.,.,treat)
	D_long = editmissing(D_long,0)
    I_long=(D_long:==0)
    st_view(X_long=.,.,(cov))
	X_long = editmissing(X_long,0)

    st_view(Y_long=.,.,outcome)
	Y_long = editmissing(Y_long,0)
	num_obs=rows(Y_long)
	p=cols(X_long)

    st_view(unit_long=.,.,unit)
	st_view(time_long=.,.,time)
	unit_long = editmissing(unit_long,0)
	time_long = editmissing(time_long,0)
	num_unit=rows(uniqrows(unit_long))
	num_time=rows(uniqrows(time_long))

    I_NT = J(num_unit,num_time,1)

    //get XX_inv (p×p dimension)
	if(p>0){
		X_long_I = X_long:*I_long //NT*p
		XX_inv = transposeonly(X_long_I)*X_long_I
		_invsym(XX_inv) //p*p
	}
	
    //reshape Y_long; unit_long; time_long
    Y_panel=rowshape(Y_long,num_unit) //N*T
    D_panel=rowshape(D_long,num_unit) //N*T
    I_panel=rowshape(I_long,num_unit) //N*T
	
    //initialization
	st_view(alpha_hat=.,.,alpha)
	alpha_hat = editmissing(alpha_hat,0)
    alpha_hat = rowshape(alpha_hat,num_unit) //N*T
	st_view(xi_hat=.,.,xi)
	xi_hat = editmissing(xi_hat,0)
    xi_hat = rowshape(xi_hat,num_unit) //N*T
	st_view(mu_hat=.,.,mu)	
	mu_hat = editmissing(mu_hat,0)
    mu_hat = rowshape(mu_hat,num_unit)  //N*T
	mc_hat = J(num_unit,num_time,0)  //N*T

    // EM Algorithm
    for(m=1;m<=iter;m++) {
        Y_hat_panel = Y_panel-mu_hat-alpha_hat-xi_hat-mc_hat
        Y_hat_panel_I = Y_hat_panel:*I_panel
        Y_hat_long_I = rowshape(Y_hat_panel_I,num_obs)
        if(p>0){
			XY = transposeonly(X_long_I)*Y_hat_long_I
			beta_hat = XX_inv*XY
		}
		
        //calculate W, a N*T matrix
        W=J(num_unit,num_time,0) //N*T
        if(p>0){
            W = (Y_panel - rowshape(X_long*beta_hat,num_unit)):*I_panel
        }else{
            W = Y_panel:*I_panel
        }
        W = W + (mu_hat  + alpha_hat + xi_hat + mc_hat):*D_panel
        
        //calculate W_mean; W_unit_mean; W_time_mean
        W_all_mean = J(num_unit,num_time,sum(W)/num_obs)
        W_unit_mean = J(num_unit,num_time,0)
        W_time_mean = J(num_unit,num_time,0)
        if(force=="two-way"|force=="unit"){ // mean across T
			W_unit_mean = rowsum(W):*I_NT/num_time
		}
		
        if(force=="two-way"|force=="time"){ // mean across N
			W_time_mean = I_NT:*colsum(W)/num_unit
		}
        if(force=="two-way"){
			W_hat = W - (W_unit_mean+W_time_mean) + W_all_mean
		}
		if(force=="unit"){
			W_hat = W - W_unit_mean
		}
		if(force=="time"){
			W_hat = W - W_time_mean
		}
		if(force=="none"){
			W_hat = W - W_all_mean
		}
		
        //svd
		if(rows(W_hat)<cols(W_hat)){
			trans=1
			W_hat = transposeonly(W_hat)
		}else{
			trans=0
		}
		W_hat = W_hat/num_obs
		svd(W_hat,U=.,svd_val=.,Vt=.)		
		svd_val=svd_val:-lambda
		svd_val=svd_val:*(svd_val:>0)
		mc_fit=U*diag(svd_val)*Vt

		if(trans==1){
			mc_fit=transposeonly(mc_fit)
		}
		
        //update
		mc_hat=mc_fit*num_obs
		if(force=="two-way"|force=="unit"){
			alpha_hat = W_unit_mean-W_all_mean
		}
		if(force=="two-way"|force=="time"){
			xi_hat = W_time_mean - W_all_mean
		}
		mu_hat = W_all_mean

        //predict
		Y_predict = mu_hat+alpha_hat+xi_hat+mc_hat
		if(p>0){
			Y_predict = Y_predict + rowshape(X_long*beta_hat,num_unit)
		}

        if(m==1){
			Y_predict_old = Y_predict
		}else{ 
			Y_predict_delta = (Y_predict-Y_predict_old):*I_panel
			crit = sqrt(sum(Y_predict_delta:*Y_predict_delta))/sqrt(sum(Y_predict:*Y_predict))
			if (crit<Tol){
				st_numscalar("r(stop)", 1)
				break
			} 
			st_numscalar("r(stop)", 0)
			Y_predict_old = Y_predict
		}
    }
	Y_predict = rowshape(Y_predict,num_obs)
    st_addvar("double", newvarname)
	st_store(., newvarname,Y_predict)

	if(p>0){
		st_matrix("r(coef)",transposeonly(beta_hat))
	}
	st_numscalar("r(cons)", sum(W_all_mean)/num_obs)

}


void val_gen(string treat, 
			 string unit, 
			 string time, 
			 string touse, 
			 real fold,
			 real period,
			 string validation,
			 real seed)
{
	real matrix panel_treat,x,panel_touse,groupid,slice_treat,slice_touse,index,targetindex
	real scalar rownum,colnum,i,j,k,groupnum,colcut,mod,num_infold
	st_view(x=.,.,(unit,time,treat,touse))

	rownum=rows(uniqrows(x[,1]))
	colnum=rows(uniqrows(x[,2]))

	panel_treat=x[,3]
	panel_treat=rowshape(panel_treat,rownum)
	panel_touse=x[,4]
	panel_touse=rowshape(panel_touse,rownum)

	groupid=J(rownum,colnum,.)
	groupnum=1

	if(period>=0.5*colnum) {
		_error("CVnobs is too large.")
	}

	if(fold<=1) {
		_error("Folds should be larger than 1.")
	}

	colcut=floor(colnum/period)

	if(colcut*period==colnum){
		mod = 1
	}
	else{
		mod = 0
	}

	for(i=1;i<=rownum;i++) {
		for(j=0;j<colcut;j++) {
			slice_treat=(panel_treat[i,(j*period+1)::((j+1)*period)]:==0)
			slice_touse=panel_touse[i,(j*period+1)::((j+1)*period)]	
			if(sum(slice_touse:*slice_treat)>0) {
				for(k=(j*period+1);k<=(j+1)*period;k++) {
					groupid[i,k]=groupnum
				}
				groupnum=groupnum+1
			}
		}
	
		if(mod==1) {
			continue
		}
	
		slice_treat=(panel_treat[i,(colcut*period+1)::colnum]:==0)
		slice_touse=panel_touse[i,(colcut*period+1)::colnum]
	
		if(sum(slice_touse:*slice_treat)>0) {
			for(k=(colcut*period+1);k<=colnum;k++) {
				groupid[i,k]=groupnum
			}
			groupnum=groupnum+1
		}
	}

	index=1::(groupnum-1)
	rseed(seed)
	_jumble(index)
	num_infold=floor((groupnum-1)/fold)

	if(num_infold==0) {
		_error("kfold is too large.")
	}

	for(i=0;i<fold;i++) {
		targetindex=index[(i*num_infold+1)::(i+1)*num_infold]
		for(j=1;j<=num_infold;j++) {
			_editvalue(groupid,targetindex[j],-(i+1))
		}
	}

	_editmissing(groupid,0)
	groupid=-groupid:*(groupid:<0)
	groupid=colshape(groupid,1)
	st_store(.,validation,groupid)
}


void eig_v(string outcome, string treat, string unit, string time)
{
	real matrix tr,use,x
	real scalar minid,maxid,rownum,num
	real matrix out_outcome,index,s

	st_view(tr=.,.,treat)
	use=(tr:==0)
	st_view(x=.,.,(unit,time,outcome))

	minid=min(x[,1])
	maxid=max(x[,1])
	rownum=maxid-minid+1
	num=rows(use)
	index=rowshape(use,rownum)
	out_outcome=editmissing(rowshape(x[,3],rownum),0):*index
	
	out_outcome=out_outcome/num
	
	s=svdsv(out_outcome)

	st_numscalar("r(max_sing)", max(s))
	
}

void p_treat(string Treat,
			 string unit, 
			 string time, 
			 string outputs,
			 real period)
{
	real matrix treat,panel_treat,x,placebotreat
	real scalar minid,maxid,mintime,maxtime,rownum,colnum,i,j,k
	st_view(treat=.,.,Treat)
	st_view(x=.,.,(unit,time))

	minid=min(x[,1])
	maxid=max(x[,1])
	mintime=min(x[,2])
	maxtime=max(x[,2])
	rownum=maxid-minid+1
	colnum=(maxtime-mintime)+1

	panel_treat=treat
	panel_treat=rowshape(panel_treat,rownum)

	placebotreat=panel_treat

	for(i=1;i<=rownum;i++){
		for(j=colnum;j>=1;j--){
			if(panel_treat[i,j]==1){
				for(k=1;k<=period;k++){
					if((j-k)>0){
						placebotreat[i,j-k]=placebotreat[i,j-k]+1
					}
				}
			}	 
		}
	}
	placebotreat=placebotreat:>=1
	placebotreat=colshape(placebotreat,1)
	st_addvar("int", outputs)
	st_store(., outputs,placebotreat)
}

/*NEW LSJ*/
/*Carryover Effect Treatment*/
void c_treat(string Treat,
			 string unit, 
			 string time, 
			 string outputs,
			 real period)
{
	real matrix treat,panel_treat,x,carryovertreat
	real scalar minid,maxid,mintime,maxtime,rownum,colnum,i,j,k
	st_view(treat=.,.,Treat)
	st_view(x=.,.,(unit,time))

	minid=min(x[,1])
	maxid=max(x[,1])
	mintime=min(x[,2])
	maxtime=max(x[,2])
	rownum=maxid-minid+1
	colnum=(maxtime-mintime)+1
	panel_treat=treat
	panel_treat=rowshape(panel_treat,rownum)

	carryovertreat=panel_treat

	for(i=1;i<=rownum;i++){
		for(j=colnum;j>=1;j--){
			if(panel_treat[i,j]==1){
				for(k=1;k<=period;k++){
					if((j+k)<=colnum){
						carryovertreat[i,j+k]=carryovertreat[i,j+k]+1
					}
				}
			}	 
		}
	}
	carryovertreat=carryovertreat:>=1
	carryovertreat=colshape(carryovertreat,1)
	st_addvar("int", outputs)
	st_store(., outputs,carryovertreat)
}
/*LSJ END*/

void waldF(string scalar Id, 
		   string scalar Time, 
		   string scalar TE, 
		   string scalar ATTs, 
		   string scalar S, 
		   real scalar from, 
		   real scalar to, 
		   string scalar to_use)
{
	real matrix panel_s,panel_atts,panel_e,id,s,time,atts,e,touse,panel_touse
	real scalar maxid,maxtime,mintime,timespan,Ot,nominator,denominator,Fobs,i,j,mi
	if(from>to){
		_error("preperiod should less than to")
	}
	
	st_view(id=.,.,Id)
	st_view(s=.,.,S)
	st_view(time=.,.,Time)
	st_view(atts=.,.,ATTs)
	st_view(e=.,.,TE)
	st_view(touse=.,.,to_use)

	maxid=max(id)
	maxtime=max(time)
	mintime=min(time)

	timespan=maxtime-mintime+1
	panel_s=s
	panel_s=rowshape(panel_s,maxid)
	panel_atts=atts
	panel_atts=rowshape(panel_atts,maxid)
	panel_e=e
	panel_e=rowshape(panel_e,maxid)
	panel_touse=touse
	panel_touse=rowshape(panel_touse,maxid)

	Ot=0
	for(i=1;i<=maxid;i++) {
		if(missing(panel_s[i,])==timespan){
			continue
		}
		mi=0
		for(j=1;j<=timespan;j++) {
			if(panel_s[i,j]>=from & panel_s[i,j]<=to & missing(panel_s[i,j])==0 & panel_touse[i,j]==1) {
				mi=mi+1
			}
		Ot=Ot+mi
		mi=0
		}
	}
	nominator=0
	denominator=0
	mi=min((abs(from)+1+to,abs(min(s))+1+to))

	for(i=1;i<=maxid;i++) {
		if(missing(panel_s[i,])==timespan) {
			continue
		}
		for(j=1;j<=timespan;j++) {
			if(panel_s[i,j]>=from & panel_s[i,j]<=to & missing(panel_s[i,j])==0 & panel_touse[i,j]==1 &missing(panel_e[i,j])==0) {
				nominator=nominator+(panel_e[i,j]^2-(panel_e[i,j]-panel_atts[i,j])^2)/mi
				denominator=denominator+((panel_e[i,j]-panel_atts[i,j]))^2/(Ot+mi)
			}
		}
	}
	Fobs=nominator/denominator
	st_numscalar("r(F)", Fobs)
}

void permute(string demissing_treat,string unit, string time, string permu_treat)
{
	real matrix treat,treat2,panel_treat,panel_treat2,x,index
	real scalar minid,maxid,rownum,i
	st_view(treat=.,.,demissing_treat)
	st_view(x=.,.,(unit,time))

	minid=min(x[,1])
	maxid=max(x[,1])
	rownum=maxid-minid+1
	

	panel_treat=treat
	panel_treat=rowshape(panel_treat,rownum)
	
	index=minid::maxid
	_jumble(index)
	panel_treat2=panel_treat
	
	for(i=1;i<=rownum;i++) {	
		panel_treat2[i,]=panel_treat[index[i],]
	}
	treat2=colshape(panel_treat2,1)
	st_addvar("int", permu_treat)
	st_store(., permu_treat,treat2)
}

end							   
						   
