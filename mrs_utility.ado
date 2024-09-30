********************************************************************************
*Reversing the reversal
********************************************************************************
*Utility to assess changes in MRS across transformations of the dependent variable.
********************************************************************************

cap program drop mrs_utility
program mrs_utility, rclass

	*=====================================
	*1. Initial setup
	*=====================================

	syntax, 								/// 
	numerator(string)						///
	denominator(string) 					///		 			
	[										///
	cost_type(integer 2)					///
	precision(integer 100)					///
	gen										///
	graph									///
	]	
	
	qui {
	
	*-------------------------------------
	*	1.1 Housekeeping
	*-------------------------------------
		
	*Store the previous model
	tempname prevmodel
	estimates store `prevmodel'
	
	*Instantiate some local names to be used later. 
	tempname tmp_mat 		// temporary matrix used in various parts of the program
	tempname tmp_mat2 		// temporary matrix used in various parts of the program	
	tempname tmp_w_mat 		// records the coefficients from regressions of hd. Should eventually be renamed "hd_result"
	tempname d_result		// records signs of coefficients in regressions of hd
	
	*-------------------------------------
	*	1.2 Get original quantities needed later
	*-------------------------------------
	
	*Get the original estimates
	tempname orig_b
	matrix `orig_b' = e(b)
	tempname orig_b_backup
	matrix `orig_b_backup' = e(b)
			
	*Get the full original command and the dependent variable
	local full_command  "`e(cmdline)'"
	local depvar "`e(depvar)'"
	
	*Get a rescaled version of the original variable
	sum `depvar', meanonly
	local min_depvar = r(min)
	local max_depvar = r(max)
	
	if "`range'" == "" {
		local scale_min = `min_depvar'
		local scale_max =  `max_depvar'		
	}
	else {
		tokenize `range'
		local scale_min = `1'
		local scale_max =  `2'
	}
		
	tempvar depvar_rescaled	
	gen `depvar_rescaled' = (((`depvar'-`min_depvar')/(`max_depvar'-`min_depvar'))*(`scale_max'-`scale_min')) + `scale_min' // division normalises to [0,1], multiplication stretches to the width of the desired scale, addition shifts to the start point of the desired scale.  
	
	*Get the labels into a matrix (needed for the python routine)
	tempname _labels_depvar_rescaled 								// can't use a temporary matrix here because python won't recognise it in that case. 
	levelsof `depvar_rescaled', matrow(_labels_depvar_rescaled)
	
	*=====================================
	*2. Assess if reversals are at all possible
	*=====================================
	
	*-------------------------------------
	*2.1 Get hd 
	*-------------------------------------
	
	levelsof `depvar', local(levels_depvar) // gets levels of depvar. Also used in later sections. 
	local n = 1 							// initialises counter (since depvar may not start at 1.)	
	foreach i of local levels_depvar {
		tempvar dstub`n'
		gen 	`dstub`n'' = 0 if `depvar' != . & e(sample)==1
		qui replace `dstub`n'' = 1 if `depvar' <= `i'	& e(sample)==1
		if "`dstub'" != "" gen `dstub'`n' =`dstub`n''
		local ++n
	}
	local n_levels_depvar = `n' -1

	*-------------------------------------
	*2.2 Run regressions of hd 
	*-------------------------------------

	local nrows_d_result = `n' - 2 	// number of regressions to be run. -2 because the local n, defined above, increments one more time than is intuitive. 
	forvalues n=1(1)`nrows_d_result' {
		
		*Get the estimation command
		local to_run = subinword("`full_command'", "`depvar'", "`dstub`n''", 1) 	// changes the dependent variable of the command originally run.
		
		*Run the regression. 
		`to_run' 
		
		*Record the result
		matrix `tmp_mat' = e(b) 
		if `n'==1 matrix `tmp_w_mat'=`tmp_mat' 
		else	  matrix `tmp_w_mat' = (`tmp_w_mat' \ `tmp_mat')
		
		*Just get the signs	and place in matrix
		mata : st_matrix("`tmp_mat2'", st_matrix("`tmp_mat'") :> 0)		
		if `n'==1 matrix `d_result'=`tmp_mat2'
		else	  matrix `d_result' = (`d_result' \ `tmp_mat2')
	
	}
		
	*=====================================
	*3. Implement the cost-function approach
	*=====================================

	*-------------------------------------
	*3.1 Some setups
	*-------------------------------------
			
	local explanatory_vars "`:colnames `orig_b''"
	local explanatory_vars = subinstr("`explanatory_vars'", "_cons", "",1)

	matselrc `tmp_w_mat' _numerator, c(`numerator')
	matselrc `tmp_w_mat' _denominator, c(`denominator')
	matrix _bd_coeffs = (_numerator, _denominator) // For getting the coefficients on hd imported into python.

	*-------------------------------------
	*3.2 Check whether signs of denom stay the same. If no, say that the ratio is bounded. 
	*-------------------------------------

	matrix colnames `d_result'=`:colnames `tmp_w_mat'' 
	tempname signs_denominator
	matselrc `d_result' `signs_denominator', c(`denominator')
	
	tempname is_reversible_mat
	matrix `is_reversible_mat' = J(1,rowsof(`signs_denominator'),1)*`signs_denominator' // computes column wise sum 
	
	*Irreversible when the columnwise sum is zero or equal to the number of coefficients:
	if `is_reversible_mat'[1,1] != 0 & `is_reversible_mat'[1,1] != `nrows_d_result' {
		noi dis as error "The denominator of the coefficient ratio is reversible. The coefficient ratio is therefore in principle unbounded and this routine will likely provide incorrect results."	
	}
		
	*-------------------------------------
	*3.3 Get costs
	*-------------------------------------
	
	cap matrix drop result_min
	cap matrix drop result_max
	
	*noi python script "/Users/casparkaiser/Library/CloudStorage/OneDrive-Personal/documents/oxford_phd/reverse_reversal_project/do/do_ck/mrs_utility_python_v250924.py" 	
	python script "`c(sysdir_plus)'py/mrs_utility_python.py"			
	
	local ratio_0 : di %4.3f `ratio_0'
	local min_ratio : di %4.3f `min_ratio'
	local max_ratio : di %4.3f `max_ratio'
	
	noi dis "No transformation: b(`numerator')/(`denominator')=`ratio_0'"
	noi dis "Min:               b(`numerator')/(`denominator')=`min_ratio'"
	noi dis "Max:               b(`numerator')/(`denominator')=`max_ratio'"
	
	if "`gen'" == "gen" {
	*For minimisation results
		cap drop _tmp1 
		cap drop _tmp2
		cap drop cost_min 
		cap drop ratio_min
		
		svmat result_min, names(_tmp)
		rename _tmp1 cost_min
		rename _tmp2 ratio_min

		*For maximisation results
		cap drop _tmp1 
		cap drop _tmp2		
		cap drop cost_max 
		cap drop ratio_max
		
		svmat result_max, names(_tmp)
		rename _tmp1 cost_max
		rename _tmp2 ratio_max
	}
	if "`graph'" == "graph" {
		if "`gen'" != "gen" noi dis as error "The graph option can only be used when the gen option is specified."
		else {
			twoway (line ratio_min cost_min) (line ratio_max cost_max), ///
			legend(order(1 "Minimum ratio" 2 "Maximum ratio") pos(6) cols(2)) ///
			ytitle("Coefficient ratio") xtitle("Cost") title("Coefficient ratio of {bf:`numerator'} to {bf:`denominator'} at given costs.") ///
			plotregion(fcolor(white)) yscale(lcolor(black)) xscale(lcolor(black)) ///
			xlabel(0(0.1)1, format(%4.1f) nogrid glcolor(black) glpattern(dot)) ///
			ylabel(, glcolor(black) glpattern(dot)) //
		}
	}
			
	*=====================================
	*5. Restore the original model.
	*=====================================
	
	matrix drop _labels_depvar_rescaled
	matrix drop _bd_coeffs
	matrix drop _denominator
	matrix drop _numerator

	qui estimates restore `prevmodel'
	
	}	// ends the qui condition
	
end
