********************************************************************************
*Reversing the reversal
********************************************************************************
*Utility to search for reversals and changes in statistical significance. 
********************************************************************************


cap program drop rr_utility
program rr_utility, rclass

	*=====================================
	*1. Initial setup
	*=====================================

	syntax, [ 								/// 
	start(real -2) end(real 2) 				/// Specifies the smallest and largest value of c over which should be searched.
	prec(real 0.1) 							/// Specifies the 'density' or 'precision' of the grid of c. For example prec(0.1) says that we evaluate values of c in steps of 0.1.
	critval(real 0.05) 						/// Specifies the alpha level where we speak of statistical significance. 
	fast 									/// Runs much faster but does not compute critical values for p-values. 
	PYthon 				 						/// Specifies that the least non-linear reversal condition is searched for. Requires that Python can be called from within Stata.  
	scale_min(real 99) scale_max(real 99) 	/// Specifies the minimum and maximum we want the scale to always be on. Not currently fully tested. May act weird. 
	GAmma(real 0) 							/// Specifies gamma
	keep(string) 							/// Specifies list of variables to be kept in the displayed results table(s).
	transpose 								/// Alternative layout for the results table. Here, rows are variables.
	dstub(string) 							/// Specifies that the binary dummy should be saved and storted in a stub specified by string. 
	]		
	
	qui {
	
	*-------------------------------------
	*	1.1 Housekeeping
	*-------------------------------------
	
	*Store the previous model
	tempname prevmodel
	estimates store `prevmodel'
	
	*Instantiate some local names to be used later. 
	tempname b_result 		// records (transformed) coefficients
	tempname p_result 		// records p-values from transformed regressions
	tempname hd_p_vals		// records p-values from regressions of hd
	tempname tmp1 			// temporary matrix used in various parts of the program
	tempname tmp2			// temporary matrix used in various parts of the program
	tempname tmp_mat 		// temporary matrix used in various parts of the program
	tempname tmp_mat2
	tempname d_result 		// records which regressions of hd are responsible for reversal
	tempname dlabels 		// record differences in labels. Needed for "fast" routine
	tempname tmp_w_mat 		// records the coefficients from regressions of hd. Should eventually be renamed "hd_result"
	tempname tmp_w_mat_pos 	// records sum of positive coefficients in regressions of hd
	tempname tmp_w_mat_neg 	// records sum of negative coefficients in regressions of hd
	tempname tmp_bin_mat 	// records coefficients of regressions of bin (bin=variable recording whether respondent answer with category k; needed for dq computations)
	tempname tmp_pr_mat 	// records predicted probability of being in a given category (needed for dq computations)
	tempname b_gamma_result	// records coefficients shifted by gamma
	tempname p_gamma_result // records p-values shifted by gamma
	
	*-------------------------------------
	*	1.2 Get original quantities needed later
	*-------------------------------------
	
	*Get the original estimates
	tempname orig_b
	matrix `orig_b' = e(b)
	tempname orig_b_backup
	matrix `orig_b_backup' = e(b)
	
	tempname orig_v
	matrix `orig_v' = e(V)
	
	*Get N and number of parameters P
	local P = colsof(e(b))
	local N = e(N)
	
	*Get the full original command and the dependent variable
	local full_command  "`e(cmdline)'"
	local depvar "`e(depvar)'"
	
	*Get the signs of the original estimates
	tempname signs_mat signs_mat_pos signs_mat_neg
	mata : st_matrix("`signs_mat_pos'", st_matrix("`orig_b'") :> 0) 		// everything that is positive is 1
	mata: st_matrix("`signs_mat_neg'", (st_matrix("`orig_b'") :< 0) :* -1)  // everything that is negative is -1
	matrix `signs_mat' = `signs_mat_pos' + `signs_mat_neg' 					// adds the two to get positive and negative signs. Working with Stata matrices is tedious. 
	
	*Get the gamma mat
	tempname gamma_mat
	matrix `gamma_mat' = `gamma' * `signs_mat'
	
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
	*2.2 Run regressions of hd and check whether signs stay the same
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
		
		*Save p-values
		tempname hd_p_vals_tmp
		tempname hd_covariance
		matrix `hd_covariance' = vecdiag(e(V))
		 
		mata : st_matrix("`hd_p_vals_tmp'", 2*(J(1,`P',1)-t((`N'-`P'), abs(st_matrix("`tmp_mat'") :/ sqrt(st_matrix("`hd_covariance'")))))) // not sure why I put a capture here. 
		if `n'==1 matrix `hd_p_vals'=`hd_p_vals_tmp' 
		else	  matrix `hd_p_vals' = (`hd_p_vals' \ `hd_p_vals_tmp')				
	}
	
	*-------------------------------------
	*2.3 Get a matrix that simply indicates reversibility (0=No; 1=Yes)
	*-------------------------------------
	
	tempname is_reversible_part1
	matrix `is_reversible_part1' = J(1,rowsof(`d_result'),1)*`d_result' // computes column wise sum and returns this to the matrix
	*Case 1: Irreversible when the columnwise sum is zero:
	tempname is_reversible_part2	
	mata : st_matrix("`is_reversible_part2'", st_matrix("`is_reversible_part1'") :== 0) 		
	*Case 2: Irreversible when the columnwise sum is == rowsof(d_result):
	tempname is_reversible_part3	
	mata : st_matrix("`is_reversible_part3'", st_matrix("`is_reversible_part1'") :== `nrows_d_result') 		// everything that is positive is 1
	
	*Now sum the two
	tempname is_reversible_part4	
	matrix `is_reversible_part4' = `is_reversible_part2' + `is_reversible_part3'

	*Make 0=1 and 1=0
	tempname is_reversible	
	mata : st_matrix("`is_reversible'", st_matrix("`is_reversible_part4'") :== 0) 			
	
	*=====================================
	*3. Implement the cost-function approach for the NON-gamma-shifted version
	*=====================================

	if "`python'"=="python" {

		*-------------------------------------
		*3.1 Some setups
		*-------------------------------------
				
		local explanatory_vars "`:colnames `orig_b''"
		local explanatory_vars = subinstr("`explanatory_vars'", "_cons", "",1)
		python clear
		
		*-------------------------------------
		*3.2 Loop over explnatory variable
		*-------------------------------------
		
		preserve
		clear 
		svmat `tmp_w_mat', names("bd")
 
		set obs `n_levels_depvar'
			
		local n = 1
		foreach explanatory_var of local explanatory_vars {
			
			*-------------------------------------
			*3.3 Prepare things for python
			*-------------------------------------

			*Skip variable if no reversal can be achieved (i.e. when all regressions of hd have the same sign.)
			cap drop neg_sign_bd
			gen neg_sign_bd = bd`n' < 0 
			replace neg_sign_bd = . if bd`n'==.
			sum neg_sign_bd, meanonly
			if r(mean)==0 | r(mean)==1  {
				gen new_labels`n' = .
				gen cost`n' = .
				local ++n
				continue
			}
			
			*Make sure that these variables, which Python will create, don't already exist. 
			cap drop python_labels
			cap drop python_cost
			
			*Get a sign variable to be imported into Python.
			cap drop sign
			gen sign = `signs_mat'[1,`n']
			
			*Get a variable with the coefficients on hd to be imported into Python.
			cap drop bd
			gen bd = bd`n'

			*-------------------------------------
			*3.4 Do the substantive Python bits
			*-------------------------------------
			
			local is_this_gamma = 0 // tells python that it should NOT actually use the value of gamma, and instead set it to zero.
			python script "`c(sysdir_plus)'py/cost_minimizer_gamma.py"			
		
			*-------------------------------------
			*3.5 Get results into the right variables
			*-------------------------------------
			
			gen new_labels`n' = python_labels
			gen cost`n' = python_cost

			*-------------------------------------
			*3.6 Iterate counter
			*-------------------------------------

			local ++n		
		}
		
		*-------------------------------------
		*3.7 Save results to matrix
		*-------------------------------------
		
		local m = `n' -1
		local costs
		forvalues j = 1(1)`m' {
			local costs `costs' cost`j'
		}
		tempname cost_mat
		mkmat `costs', matrix(`cost_mat')
		matrix `cost_mat' = `cost_mat'[1,1...]
		
		local labels
		forvalues j = 1(1)`m' {
			local labels `labels' new_labels`j'
		}
		tempname label_mat
		mkmat `labels', matrix(`label_mat')
		
		restore 
	}


	*=====================================
	*4. Implement the cost-function approach for the gamma-shifted version
	*=====================================

	if "`python'"=="python" & `gamma'!=0 {
		
		*-------------------------------------
		*4.1 Some setups
		*-------------------------------------
		
		*Get matrix with only the positive coefficients. Then find columnwise maximum
		tempname tmp_w_mat_pos 
		tempname tmp_w_mat_pos_max
		tempname tmp_w_mat_pos_min		
		mata: st_matrix("`tmp_w_mat_pos'", (st_matrix("`tmp_w_mat'") :> 0) :* st_matrix("`tmp_w_mat'"))	
		mata: st_matrix("`tmp_w_mat_pos_max'", colmax(st_matrix("`tmp_w_mat_pos'")))	
		mata: st_matrix("`tmp_w_mat_pos_min'", colmin(st_matrix("`tmp_w_mat_pos'")))	
		
		*Get matrix with only the negative coefficients. Then find columnwise minimum 
		tempname tmp_w_mat_neg
		tempname tmp_w_mat_neg_max		
		tempname tmp_w_mat_neg_min
		
		mata: st_matrix("`tmp_w_mat_neg'", (st_matrix("`tmp_w_mat'") :< 0) :* st_matrix("`tmp_w_mat'"))		
		mata: st_matrix("`tmp_w_mat_neg_max'", colmax(st_matrix("`tmp_w_mat_neg'")))	
		mata: st_matrix("`tmp_w_mat_neg_min'", colmin(st_matrix("`tmp_w_mat_neg'")))	
	
		local explanatory_vars "`:colnames `orig_b''"
		local explanatory_vars = subinstr("`explanatory_vars'", "_cons", "",1)
		python clear
		
		*-------------------------------------
		*4.2 Loop over explnatory variable
		*-------------------------------------
		
		preserve
		clear 
		svmat `tmp_w_mat', names("bd")
		set obs `n_levels_depvar'
			
		local n = 1
		foreach explanatory_var of local explanatory_vars {
			
			*-------------------------------------
			*4.3 Prepare things for python
			*-------------------------------------
			
			**If cost is zero anyway, we just set it without any computation
			*Original coefficient is positive
			if `orig_b'[1,`n'] > 0 {
				if `orig_b'[1,`n'] + `gamma' < 0 {
					gen new_labels`n' = .
					gen cost`n' = 0
					local ++n
					continue
				}
			}
			*Original coefficient is negative
			if `orig_b'[1,`n'] < 0 {
				if `orig_b'[1,`n'] - `gamma' > 0 {
					gen new_labels`n' = .
					gen cost`n' = 0
					local ++n
					continue
				}
			}
						
			**Skip variable if no reversal can be achieved even after accounting for a gamma shift that makes reversals easier.  
			*First we check if all coefficients on bd are the same. Only then will we have to do a further check.
			cap drop neg_sign_bd
			gen neg_sign_bd = bd`n' < 0 
			replace neg_sign_bd = . if bd`n'==.
			sum neg_sign_bd, meanonly
			
			*Case 1: Original b is positive
			if r(mean)==1 { // Ok, now we know that the original coefficient was (1) positive and (2) irreversible without a gamma shift. Now it's time to do a further check. 
				if abs((`n_levels_depvar'-1)*`tmp_w_mat_neg_max'[1,`n']) > abs(`gamma')  { 					
					gen new_labels`n' = .
					gen cost`n' = .
					local ++n
					continue
				}
			}
			*Case 2: Original b is negative
			if r(mean)==0 { // Ok, now we know that the original coefficient was (1) negative and (2) irreversible without a gamma shift. Now it's time to do a further check. 
				if abs((`n_levels_depvar'-1)*`tmp_w_mat_pos_min'[1,`n']) > abs(`gamma') {  					
					gen new_labels`n' = .
					gen cost`n' = .
					local ++n
					continue
				}
			}
			
			**Skip variable if no reversal can be achieved due to a gamma shift that makes reversals harder
			*Case 1: Original b is positive
			if `orig_b'[1,`n'] > 0 {
				if (`n_levels_depvar'-1)*`tmp_w_mat_pos_max'[1,`n'] < `gamma' {
					gen new_labels`n' = .
					gen cost`n' = .
					local ++n
					continue					
				}	
			}
			*Case 2: Original b is negative
			if `orig_b'[1,`n'] < 0 {
				if (`n_levels_depvar'-1)*`tmp_w_mat_neg_min'[1,`n'] > -`gamma' {
					gen new_labels`n' = .
					gen cost`n' = .
					local ++n
					continue					
				}	
			}			

									
			*Make sure that these variables, which Python will create, don't already exist. 
			cap drop python_labels
			cap drop python_cost
			
			*Get a sign variable to be imported into Python.
			cap drop sign
			gen sign = `signs_mat'[1,`n']
			
			*Get a variable with the coefficients on hd to be imported into Python.
			cap drop bd
			gen bd = bd`n'

			*-------------------------------------
			*4.4 Do the substantive Python bits
			*-------------------------------------

			local is_this_gamma = 1	// tells python that it should actually use the value of gamma		
			python script "`c(sysdir_plus)'py/cost_minimizer_gamma.py"			
		
			*-------------------------------------
			*4.5 Get results into the right variables
			*-------------------------------------
			
			gen new_labels`n' = python_labels
			gen cost`n' = python_cost
			replace cost`n' = . if cost`n' < 0.000000001 // Handles the case where base categories are not caught earlier
			
			*-------------------------------------
			*4.6 Iterate counter
			*-------------------------------------

			local ++n		
		}
		
		*-------------------------------------
		*4.7 Save results to matrix
		*-------------------------------------
		
		local m = `n' -1
		local costs
		forvalues j = 1(1)`m' {
			local costs `costs' cost`j'
		}
		tempname cost_mat_gamma
		mkmat `costs', matrix(`cost_mat_gamma')
		matrix `cost_mat_gamma' = `cost_mat_gamma'[1,1...]
		
		local labels
		forvalues j = 1(1)`m' {
			local labels `labels' new_labels`j'
		}
		tempname label_mat_gamma
		mkmat `labels', matrix(`label_mat_gamma')
		
		restore 
	}
	
	*=====================================
	*5. Loop over levels of c to find reversal conditions
	*=====================================

	*-------------------------------------
	*5.0 Get a rescaled version of the original variable
	*-------------------------------------
	
	tempvar depvar_rescaled
	sum `depvar', meanonly
	local min_depvar = r(min)
	local max_depvar = r(max)
	
	if `scale_min'==99 local scale_min = r(min)
	if `scale_max'==99 local scale_max = r(max)
	
	gen `depvar_rescaled' = (((`depvar'-`min_depvar')/(`max_depvar'-`min_depvar'))*(`scale_max'-`scale_min')) + `scale_min' // division normalises to [0,1], multiplication stretches to the width of the desired scale, addition shifts to the start point of the desired scale.  
		
	*-------------------------------------
	*5.1 Initialising things 
	*-------------------------------------

	*Tell the user how many regressions need to be run
	local N_c = ((`end' - `start')/`prec')+1
	noisily dis ""
	noisily dis "Total number of c-values to check: `N_c'. Progress:"
	
	*A counter to keep track of things
	local iter = 1
	
	*-------------------------------------
	*5.2 Begin the loop 
	*-------------------------------------
	
	forvalues c=`start'(`prec')`end'{
		
		*-------------------------------------
		*5.3 Get minus sign before the exp() function for negative c
		*-------------------------------------

		if `c' < 0 	local minus "-"
		if `c' > 0	local minus ""
		
		*=====================================
		*5. Run the fast routine
		*=====================================

		if "`fast'"=="fast" {
			
			*-------------------------------------
			*6.0 Get a rescaled version of the original variable
			*-------------------------------------

			tempvar trans_depvar_rescaled
			if (`c' < 0.0000001 & `c' > -0.0000001) gen `trans_depvar_rescaled' = `depvar_rescaled' 
			else {
				gen `trans_depvar_rescaled' = ((`minus'exp(`depvar'*`c') - `minus'exp(`min_depvar'*`c')) / (`minus'exp(`max_depvar'*`c') - `minus'exp(`min_depvar'*`c')))*(`scale_max'-`scale_min') + `scale_min'
			}

			*-------------------------------------
			*6.1 Get a matrix to record differences in labels
			*-------------------------------------

			matrix `dlabels' = J(`nrows_d_result',`nrows_d_result',0)
			local n = 0
			levelsof `trans_depvar_rescaled', local(levels_transdepvar_rescales)
			foreach l of local levels_transdepvar_rescales {
				if `n'==0 {
					local k = `l' 	// k here denotes the previous label. Do nothing in the first iteration and just update k for the next iteration.
				}
				else { 				// produces a diagonal matrix with the right dimensions to multiply.
					capture matrix `dlabels'[`n',`n'] = `k'-`l'  									// Again weird condition because of rounding errors from floating points.
					local k = `l' 	// Store l in k for the next iteration. 
				}
				local ++n
			}
			
			*-------------------------------------
			*6.2 Use the results from regressions of hd to arrive at the transformed coefficients. 
			*-------------------------------------

			*Take the dmat result and multiply each element by the matrix recording differences of labels
			matrix `tmp1' = `dlabels' * `tmp_w_mat'
			
			*Take the column-wise sum:
			matrix `tmp1' = J(1,rowsof(`tmp1'),1)*`tmp1' 
			
			*Store result in the b_result matrix
			if `c'==`start' matrix `b_result' = (`tmp1')
			if `c'!=`start' matrix `b_result' = (`b_result' \ `tmp1')
			
			*Store results when shifted by gamma in the b_result matrix
			tempname tmp_b_gamma
			matrix `tmp_b_gamma' = `tmp1' + `gamma_mat'
			if `c'==`start' matrix `b_gamma_result' = (`tmp_b_gamma') 
			if `c'!=`start' matrix `b_gamma_result' = (`b_gamma_result' \ `tmp_b_gamma')
		}
						
		*=====================================
		*7. Run the slow routine (only needed when interest is in p-values.)
		*=====================================

		if "`fast'" == "" & "`dq'"=="" {
			
			*-------------------------------------
			*7.1 Create transformed depvar and run the regression
			*-------------------------------------
			
			*Create transformed depvar			
			tempvar trans_depvar_rescaled
			if (`c' < 0.0000001 & `c' > -0.0000001) gen `trans_depvar_rescaled' = `depvar_rescaled' 
			else {
				gen `trans_depvar_rescaled' = ((`minus'exp(`depvar'*`c') - `minus'exp(`min_depvar'*`c')) / (`minus'exp(`max_depvar'*`c') - `minus'exp(`min_depvar'*`c')))*(`scale_max'-`scale_min') + `scale_min'
			}
			*Get the command to run
			local to_run = subinword("`full_command'", "`depvar'", "`trans_depvar_rescaled'", 1) 	// changes the dependent variable of the command originally run.
			
			*Run the regression
			`to_run'
			
			*Save results
			matrix `tmp1' = e(b)
			matrix `tmp2' = vecdiag(e(V))	
						
			*-------------------------------------
			*7.2 Get p-values and put coefficients into matrices  
			*-------------------------------------
			
			*Put coefficients into a results-matrix		
			if `c'==`start' matrix `b_result' = (`tmp1')			
			if `c'!=`start' matrix `b_result' = (`b_result' \ `tmp1')

			*Get p-vals and put into a results-matrix
			tempname tmp_p
			capture mata : st_matrix("`tmp_p'", 2*(J(1,`P',1)-t((`N'-`P'), abs(st_matrix("`tmp1'") :/ sqrt(st_matrix("`tmp2'")))))) // not sure why I put a capture here. 

			if `c'!=`start' matrix `p_result' = (`p_result' \ `tmp_p')
			if `c'==`start' matrix `p_result' = (`tmp_p')

			*-------------------------------------
			*7.3 Get gamma-shifted p-values and put gamma-shifted coefficients into matrices  
			*-------------------------------------
			
			*Put gamma shifted coefficients into a results-matrix		
			tempname tmp_b_gamma
			matrix `tmp_b_gamma' = `tmp1' + `gamma_mat'
			if `c'==`start' matrix `b_gamma_result' = (`tmp_b_gamma') 
			if `c'!=`start' matrix `b_gamma_result' = (`b_gamma_result' \ `tmp_b_gamma')
			
			*Get gamma shifted p-vals and put into a results-matrix
			tempname tmp_p_gamma
			capture mata : st_matrix("`tmp_p_gamma'", 2*(J(1,`P',1)-t((`N'-`P'), abs(st_matrix("`tmp_b_gamma'") :/ sqrt(st_matrix("`tmp2'")))))) // not sure why I put a capture here. 

			if `c'!=`start' matrix `p_gamma_result' = (`p_gamma_result' \ `tmp_p_gamma')
			if `c'==`start' matrix `p_gamma_result' = (`tmp_p_gamma')
			
		}
		
		*Display progress
		noisily _dots `iter' 0
		local ++iter
	}
		
	*=====================================
	*8. Use transformed estimates from any of the routines to find points of reversal across c
	*=====================================
	
	preserve
	clear
	
	svmat double `b_result', names(variable)
	qui gen tmpvar = _n
	qui tsset tmpvar
	foreach var of varlist variable1-variable`P' {
		if "`var'" == "tmpvar" continue, break  
		qui replace `var' = `var'>0 		// 1 if positive, 0 otherwise
		qui gen double new`var' = d.`var' 	// missing in the first instance, then -1 if pos to neg, and 1 if neg to pos, 0 if the same.
	}
	keep newvariable1-newvariable`P'
	tempname r_result
	mkmat newvariable1-newvariable`P', matrix(`r_result')
	
	*=====================================
	*9. Use p-value estimates from to find points of significant reversal across c (only for slow routine)
	*=====================================
	
	if "`fast'"=="" {
		clear
		svmat double `p_result', names(variable)	
		qui gen tmpvar = _n
		qui tsset tmpvar	
		foreach var of varlist variable1-variable`P' {
			if "`var'" == "tmpvar" continue, break  
			qui replace `var' = (`var'-`critval')<0 	// 1 if significant, 0 if insifnificant
			qui gen double new`var' = d.`var' 			// missing in the first instance, then -1 if significant to insignificant, and 1 if insignificant to significant, 0 if the same.
		}		
		keep newvariable1-newvariable`P'
		tempname y_result
		mkmat newvariable1-newvariable`P', matrix(`y_result')
	}
	restore
	
	*=====================================
	*10. Use transformed gamma shifted estimates from to find points of reversal across c
	*=====================================

	*-------------------------------------
	*10.1 Expand the original sign mat
	*-------------------------------------

	tempname signs_mat_expanded
	local rows_to_do = rowsof(`b_gamma_result')
	forvalues i=1(1)`rows_to_do' {
		if `i' == 1 matrix `signs_mat_expanded' = `signs_mat'
		else 		matrix `signs_mat_expanded' = `signs_mat_expanded' \ `signs_mat'
	}
		
	*-------------------------------------
	*10.2 Get a signs mat for the transformed coefficients
	*-------------------------------------

	*Get the signs of the original estimates
	tempname signs_mat_trans_gamma signs_mat_trans_gamma_pos signs_mat_trans_gamma_neg
	mata : st_matrix("`signs_mat_trans_gamma_pos'", st_matrix("`b_gamma_result'") :> 0) 		// everything that is positive is 1
	mata: st_matrix("`signs_mat_trans_gamma_neg'", (st_matrix("`b_gamma_result'") :< 0) :* -1)  // everything that is negative is -1
	matrix `signs_mat_trans_gamma' = `signs_mat_trans_gamma_pos' + `signs_mat_trans_gamma_neg' 					// adds the two to get positive and negative signs. Working with Stata matrices is tedious. 
	
	*-------------------------------------
	*10.3 Subtract the two matrices
	*-------------------------------------
	
	tempname reversal_mat_gamma
	matrix `reversal_mat_gamma' = `signs_mat_expanded' - `signs_mat_trans_gamma'

	*-------------------------------------
	*10.4 Find row with non-zero element that has smallest c
	*-------------------------------------

	/*
	Gameplan:
	Step 0: Attach the c value to the matrix 
	Step 1: Only keep the matrix for positive c. Find first non-zero element.
	Step 2: Only keep the matrix for negative c. Take absolute value of c, sort on c. Find first 1.
	Step 3: Combine and find smallest c. 
	*/
	
	*Step 0: 
	local n = 1
	tempname c_vals
	matrix `c_vals' = J(rowsof(`reversal_mat_gamma'),1,.)
	forvalues c=`start'(`prec')`end'{
		capture matrix `c_vals'[`n',1] = `c' // not sure why there's a capture here. 
		local ++n
	}
	
	matrix `reversal_mat_gamma' = (`reversal_mat_gamma', `c_vals')
	local c_num = colsof(`reversal_mat_gamma') 
	local c_num1 = `c_num' - 1
	local c_num2 = `c_num' - 2	
	
	*Step 1:
	preserve
	clear
	svmat `reversal_mat_gamma', names(tmp)
	keep if tmp`c_num' >= 0

	forvalues j=1/`c_num2' {
		sum tmp`c_num' if tmp`j'!=0
		replace tmp`j' = r(min) 	
	}
	keep if _n==1
	tempname min_cb_gamma_value1
	mkmat tmp1-tmp`c_num2', matrix(`min_cb_gamma_value1')
	restore
	
	
	*Step 2:
	preserve
	clear
	svmat `reversal_mat_gamma', names(tmp)
	keep if tmp`c_num' < 0
	replace tmp`c_num' = abs(tmp`c_num')
	sort tmp`c_num'
	forvalues j=1/`c_num2' {
		sum tmp`c_num' if tmp`j'!=0
		replace tmp`j' = r(min) 	
	}
	keep if _n==1
	tempname min_cb_gamma_value2
	mkmat tmp1-tmp`c_num2', matrix(`min_cb_gamma_value2')
	restore

	*Step 3:
	tempname combined
	matrix `combined' = `min_cb_gamma_value1' \ `min_cb_gamma_value2'
	preserve 
	clear
	svmat `combined', names(tmp)
	gen sign = 1 in 1
	replace sign = -1 in 2
	
	forvalues j=1/`c_num2' {
		sum tmp`j' 
		sum sign if tmp`j' == r(min)
		local sign = r(mean)
		sum tmp`j' 		
		replace tmp`j' = r(min) * `sign'	
	}
	keep if _n==1
	tempname min_cb_gamma
	mkmat tmp1-tmp`c_num2', matrix(`min_cb_gamma_value2')
	restore	
	
	*-------------------------------------
	*10.5 OLD ROUTINE to get the r(rgamma) matrix
	*-------------------------------------
		
	preserve
	clear
	
	svmat double `b_gamma_result', names(variable)
	qui gen tmpvar = _n
	qui tsset tmpvar
	foreach var of varlist variable1-variable`P' {
		if "`var'" == "tmpvar" continue, break  
		qui replace `var' = `var'>0 		// 1 if positive, 0 otherwise
		qui gen double new`var' = d.`var' 	// missing in the first instance, then -1 if pos to neg, and 1 if neg to pos, 0 if the same.
	}
	keep newvariable1-newvariable`P'
	tempname r_gamma_result
	mkmat newvariable1-newvariable`P', matrix(`r_gamma_result')
	
	*=====================================
	*11. Use gamma shifted p-value estimates from to find points of significant reversal across c (only for slow routine)
	*=====================================
	
	if "`fast'"=="" {
		clear
		svmat double `p_gamma_result', names(variable)	
		qui gen tmpvar = _n
		qui tsset tmpvar	
		foreach var of varlist variable1-variable`P' {
			if "`var'" == "tmpvar" continue, break  
			qui replace `var' = (`var'-`critval')<0 	// 1 if significant, 0 if insifnificant
			qui gen double new`var' = d.`var' 			// missing in the first instance, then -1 if significant to insignificant, and 1 if insignificant to significant, 0 if the same.
		}		
		keep newvariable1-newvariable`P'
		tempname y_gamma_result
		mkmat newvariable1-newvariable`P', matrix(`y_gamma_result')
	}
	restore	

	*=====================================
	*12. Clean-up matrices
	*=====================================

	*-------------------------------------
	*12.1 Get the original variable names as column names 
	*-------------------------------------

	if "`python'"!="" matrix colnames `label_mat'=`explanatory_vars'
	if "`python'"!="" matrix colnames `cost_mat'=`explanatory_vars'
		
	matrix colnames `b_result'=`:colnames `tmp_w_mat'' 
	matrix colnames `r_result'=`:colnames `tmp_w_mat'' 
	matrix colnames `r_gamma_result'=`:colnames `tmp_w_mat'' 
	
	matrix colnames `d_result'=`:colnames `tmp_w_mat'' 
	matrix colnames `hd_p_vals'=`:colnames `tmp_w_mat'' 	
	if "`fast'"=="" matrix colnames `p_result'=`:colnames `tmp_w_mat'' 
	if "`fast'"=="" matrix colnames `y_result'=`:colnames `tmp_w_mat'' 
	if "`fast'"=="" matrix colnames `y_gamma_result'=`:colnames `tmp_w_mat'' 

	matrix colnames `b_gamma_result'=`:colnames `tmp_w_mat'' 
	if "`fast'"==""  matrix colnames `p_gamma_result'=`:colnames `tmp_w_mat'' 
	
	*-------------------------------------
	*12.2 Attach the corresponding value of c to the results matrices
	*-------------------------------------
	
	local n = 1
	tempname c_vals
	matrix `c_vals' = J(rowsof(`b_result'),1,.)
	forvalues c=`start'(`prec')`end'{
		capture matrix `c_vals'[`n',1] = `c' // not sure why there's a capture here. 
		local ++n
	}
	
	matrix `b_result' = (`b_result', `c_vals')
	matrix `r_result' = (`r_result', `c_vals') 	
	matrix `r_gamma_result' = (`r_gamma_result', `c_vals') 	
	
	if "`fast'"=="" matrix `p_result' = (`p_result', `c_vals') 
	if "`fast'"=="" matrix `y_result' = (`y_result', `c_vals') 	
	if "`fast'"=="" matrix `y_gamma_result' = (`y_gamma_result', `c_vals') 	
	
	matrix `b_gamma_result' = (`b_gamma_result', `c_vals') 
	if "`fast'"=="" matrix `p_gamma_result' = (`p_gamma_result', `c_vals') 

	*=====================================
	*13. Preparare for main display to user
	*=====================================

	*-------------------------------------
	*13.0 Get minimum and maximum p-value
	*-------------------------------------
	
	**Get maximum p-value
	preserve
	clear
	svmat `hd_p_vals', names("tmp")
	ds
	local vars = "`r(varlist)'"
	local n = 1
	foreach var in `vars' {
		sum `var'
		replace `var' = r(max)
		replace `var' = 1 if `is_reversible'[1,`n'] == 1
		local ++n
	}
	keep if _n ==1
	tempname max_pval
	mkmat `vars', matrix(`max_pval')
	restore

	**Get minimum p-value
	preserve
	clear
	svmat `hd_p_vals', names("tmp")
	ds
	local vars = "`r(varlist)'"
	foreach var in `vars' {
		sum `var'
		replace `var' = r(min) 	
	}
	keep if _n ==1
	tempname min_pval
	mkmat `vars', matrix(`min_pval')
	restore
	
	*-------------------------------------
	*13.1 Get minimum c-value
	*-------------------------------------

	preserve
	clear
	local c_num = colsof(`r_result')
	local c_num1 = `c_num' - 1
	local c_num2 = `c_num' - 2
	
	svmat `r_result', names(tmp)
	gen orig_c = tmp`c_num' 
	
	forvalues j=1/`c_num' {
		replace tmp`j' = abs(tmp`j') 	
	}
	forvalues j=1/`c_num1' {
		sum tmp`c_num' if tmp`j'==1
		sum orig_c if tmp`c_num'==r(min) & tmp`j'==1
		replace tmp`j' = r(mean) 	
	}
	keep if _n==1
	tempname min_c_value
	mkmat tmp1-tmp`c_num1', matrix(`min_c_value')
	restore
	
	if "`fast'"=="" {

		*-------------------------------------
		*13.2 Get minimum c-value to render significant coefficient insignificant
		*-------------------------------------

		/*
		Gameplan:
		Step 1: Only keep the matrix for positive c. Find first -1.
		Step 2: Only keep the matrix for negative c. Take absolute value of c, sort on c. Find first 1.
		Step 3: Combine and find smallest c. 
		*/ 
		
		*Step 1:
		preserve
		clear
		svmat `y_result', names(tmp)
		keep if tmp`c_num' >= 0

		forvalues j=1/`c_num2' {
			sum tmp`c_num' if tmp`j'==-1
			replace tmp`j' = r(min) 	
		}
		keep if _n==1
		tempname min_cp_value1
		mkmat tmp1-tmp`c_num2', matrix(`min_cp_value1')
		restore
		
		
		*Step 2:
		preserve
		clear
		svmat `y_result', names(tmp)
		keep if tmp`c_num' < 0
		replace tmp`c_num' = abs(tmp`c_num')
		sort tmp`c_num'
		forvalues j=1/`c_num2' {
			sum tmp`c_num' if tmp`j'==1
			replace tmp`j' = r(min) 	
		}
		keep if _n==1
		tempname min_cp_value2
		mkmat tmp1-tmp`c_num2', matrix(`min_cp_value2')
		restore

		*Step 3:
		tempname combined
		matrix `combined' = `min_cp_value1' \ `min_cp_value2'
		preserve 
		clear
		svmat `combined', names(tmp)
		gen sign = 1 in 1
		replace sign = -1 in 2
		
		forvalues j=1/`c_num2' {
			sum tmp`j' 
			sum sign if tmp`j' == r(min)
			local sign = r(mean)
			sum tmp`j' 		
			replace tmp`j' = r(min) * `sign'	
		}
		keep if _n==1
		tempname min_cp_value
		mkmat tmp1-tmp`c_num2', matrix(`min_cp_value')
		restore
		
		*-------------------------------------
		*13.3 Get minimum c-value to render insignificant coefficient significant
		*-------------------------------------
		
		/*
		Gameplan:
		Step 1: Only keep the matrix for positive c. Find first 1.
		Step 2: Only keep the matrix for negative c. Take absolute value of c, sort on c. Find first -1.
		Step 3: Combine and find smallest c. 
		*/ 

		*Step 1:
		preserve
		clear
		svmat `y_result', names(tmp)
		keep if tmp`c_num' >= 0

		forvalues j=1/`c_num2' {
			sum tmp`c_num' if tmp`j'==1
			replace tmp`j' = r(min) 	
		}
		keep if _n==1
		tempname min_cp2_value1
		mkmat tmp1-tmp`c_num2', matrix(`min_cp2_value1')
		restore
		
		
		*Step 2:
		preserve
		clear
		svmat `y_result', names(tmp)
		keep if tmp`c_num' < 0
		replace tmp`c_num' = abs(tmp`c_num')
		sort tmp`c_num'
		forvalues j=1/`c_num2' {
			sum tmp`c_num' if tmp`j'==-1
			replace tmp`j' = r(min) 	
		}
		keep if _n==1
		tempname min_cp2_value2
		mkmat tmp1-tmp`c_num2', matrix(`min_cp2_value2')
		restore

		*Step 3:
		tempname combined
		matrix `combined' = `min_cp2_value1' \ `min_cp2_value2'
		preserve 
		clear
		svmat `combined', names(tmp)
		gen sign = 1 in 1
		replace sign = -1 in 2
		
		forvalues j=1/`c_num2' {
			sum tmp`j' 
			sum sign if tmp`j' == r(min)
			local sign = r(mean)
			sum tmp`j' 		
			replace tmp`j' = r(min) * `sign'	
		}
		keep if _n==1
		tempname min_cp2_value
		mkmat tmp1-tmp`c_num2', matrix(`min_cp2_value')
		restore	
	}
	*-------------------------------------
	*13.4 Get matrices into the right format, find original p-values, construct display matrix, implement the "keep()" option.
	*-------------------------------------
	
	**Get the right column number (removing the constant)		
	local npcols = colsof(`max_pval') - 1 
	
	**Cut to right number of columns
	*Minimum and maximum p-value
	matrix `max_pval' = `max_pval'[1...,1..`npcols'] 
	matrix `min_pval' = `min_pval'[1...,1..`npcols'] 
	*Minimum c-value.
	matrix `min_c_value' = `min_c_value'[1...,1..`npcols']
	*Original coefficients.
	tempname orig_b_display
	matrix `orig_b_display' = `orig_b'[1...,1..`npcols'] 
	*Original p-value (need to recompute them first)
	tempname orig_ses
	matrix `orig_ses' = vecdiag(`orig_v')
	tempname orig_p_display
	mata : st_matrix("`orig_p_display'", 2*(J(1,`P',1)-t((`N'-`P'), abs(st_matrix("`orig_b'") :/ sqrt(st_matrix("`orig_ses'")))))) // not sure why I put a capture here. 
	matrix `orig_p_display' = `orig_p_display'[1...,1..`npcols'] 

	**Assemble the matrix to be displayed
	tempname _to_display
	
	if "`python'" == "python" {
		if "`fast'"=="" {
			matrix `_to_display' = `orig_b_display' \ `cost_mat' \ `min_c_value' \ `orig_p_display' \ `max_pval' \ `min_pval' \ `min_cp_value' \ `min_cp2_value'
			matrix rownames `_to_display' = "1" "2" "3" "4" "5" "6" "7" "8"

		}
		else {
			matrix `_to_display' = `orig_b_display' \ `cost_mat' \ `min_c_value' \ `orig_p_display' \ `max_pval' \ `min_pval'	
			matrix rownames `_to_display' = "1" "2" "3" "4" "5" "6" 		
		}
	}
	else {
		if "`fast'"=="" {
			matrix `_to_display' = `orig_b_display' \ `min_c_value' \ `orig_p_display' \ `max_pval' \ `min_pval' \ `min_cp_value' \ `min_cp2_value'
			matrix rownames `_to_display' = "1" "2" "3" "4" "5" "6" "7" 

		}
		else {
			matrix `_to_display' = `orig_b_display' \ `min_c_value' \ `orig_p_display' \ `max_pval' \ `min_pval'	
			matrix rownames `_to_display' = "1" "2" "3" "4" "5"  		
		}
	}
	**Implement keep list (requires matselrc):
	if "`keep'" != "" matselrc `_to_display' `_to_display', c(`keep') 

	*=====================================
	*14. Preparare for additional gamma display to user
	*=====================================
	
	if `gamma' != 0 {
		
		*-------------------------------------
		*14.0 Original coefficient
		*-------------------------------------

		tempname orig_b_display_gamma
		matrix `orig_b_display_gamma' = `orig_b' + `gamma_mat'
		matrix `orig_b_display_gamma' = `orig_b_display_gamma'[1...,1..`npcols'] 

		*-------------------------------------
		*14.1 Assemble matrix to be displayed
		*-------------------------------------

		tempname _to_display_gamma
		if "`python'" == "python" {
			matrix `_to_display_gamma' = `orig_b_display_gamma' \ `cost_mat_gamma' \ `min_cb_gamma_value2'
			matrix colnames `_to_display_gamma'=`explanatory_vars'
			matrix rownames `_to_display_gamma' = "1" "2" "3"
		}
		else {
			matrix `_to_display_gamma' = `orig_b_display_gamma' \ `min_cb_gamma_value2'
			matrix colnames `_to_display_gamma'=`explanatory_vars'
			matrix rownames `_to_display_gamma' = "1" "2" 
		}		
		
		**Implement keep list (requires matselrc):
		if "`keep'" != "" matselrc `_to_display_gamma' `_to_display_gamma', c(`keep') 
		
	}
	
	*=====================================
	*15. Display to user if python is specified
	*====================================
	
	if "`python'" == "python" {
		*-------------------------------------
		*15.1 Main display
		*-------------------------------------
		
		if "`fast'"=="" {
			noi {
				dis ""
				dis ""		
				dis "{bf:Results:}"
				if "`transpose'" != "" { // sets alternative layout suggested by Anthony
					esttab matrix(`_to_display', fmt(3) transpose), mtitles("") ///
					collabels( ///
					"Original coefficient" ///
					"Minimum cost" ///
					"Minimum c-value"  ///
					"Original p-value" ///
					"Maximum p-value" ///
					"Minimum p-value" ///
					"Min. c-value for insig." ///
					"Min. c-value for sig." ) ///
					modelwidth(25) note("Missing values imply one of the following: `=char(13)'`=char(10)'" ///
					"(1) That no original coefficient was estimated. `=char(13)'`=char(10)'" ///
					"(2) That no reversal is possible. `=char(13)'`=char(10)'" ///
					"(3) That the minimum reversing c-value falls outside the search range. `=char(13)'`=char(10)'" ///
					"(4) That coefficients cannot be made significant or insignificant within the search range.") 
				}
				else { // standard layout
					esttab matrix(`_to_display', fmt(3)), mtitles("") ///
					varlabels( ///
					1 "Original coefficient:" ///
					2 "Minimum cost:" ///
					3 "Smallest c-value for reversal:"  ///
					4 "Original p-value:" ///
					5 "Maximum p-value:" ///
					6 "Minimum p-value:" ///
					7 "Smallest c-value for insig. coeff.:" ///
					8 "Smallest c-value for sig. coeff.:" ) ///
					varwidth(35) note("Missing values imply one of the following: `=char(13)'`=char(10)'" ///
					"(1) That no original coefficient was estimated. `=char(13)'`=char(10)'" ///
					"(2) That no reversal is possible. `=char(13)'`=char(10)'" ///
					"(3) That the minimum reversing c-value falls outside the search range. `=char(13)'`=char(10)'" ///
					"(4) That coefficients cannot be made significant or insignificant within the search range.") 
				}
			}
		}
		else {
			noi {
				dis ""
				dis ""		
				dis "{bf:Results:}"
				if "`transpose'" != "" { // sets alternative layout suggested by Anthony
					esttab matrix(`_to_display', fmt(3) transpose), mtitles("") ///
					collabels( ///
					"Original coefficient" ///
					"Minimum cost" ///
					"Minimum c-value"  ///
					"Original p-value" ///
					"Maximum p-value" ///
					"Minimum p-value") ///
					modelwidth(25) note("Missing values imply one of the following: `=char(13)'`=char(10)'" ///
					"(1) That no original coefficient was estimated. `=char(13)'`=char(10)'" ///
					"(2) That no reversal is possible. `=char(13)'`=char(10)'" ///
					"(3) That the minimum reversing c-value falls outside the search range. `=char(13)'`=char(10)'")
				}
				else { // standard layout
					esttab matrix(`_to_display', fmt(3)), mtitles("") ///
					varlabels( ///
					1 "Original coefficient:" ///
					2 "Minimum cost:" ///
					3 "Smallest c-value for reversal:"  ///
					4 "Original p-value:" ///
					5 "Maximum p-value:" ///
					6 "Minimum p-value:") ///
					varwidth(35) note("Missing values imply one of the following: `=char(13)'`=char(10)'" ///
					"(1) That no original coefficient was estimated. `=char(13)'`=char(10)'" ///
					"(2) That no reversal is possible. `=char(13)'`=char(10)'" ///
					"(3) That the minimum reversing c-value falls outside the search range. `=char(13)'`=char(10)'" ///
					"(4) That coefficients cannot be made significant or insignificant within the search range.") 
				}
			}		
		}
		*-------------------------------------
		*15.2 Additional gamma display
		*-------------------------------------
		
		if `gamma' != 0 {
			noi {
				dis ""
				dis ""		
				dis "{bf:Further results after applying shift by gamma:}"
				if "`transpose'" != "" { // sets alternative layout suggested by Anthony
					esttab matrix(`_to_display_gamma', fmt(3) transpose), mtitles("") ///
					collabels( ///
					"Untrans. gamma-shifted coeff." ///
					"Minimum cost" ///
					"Minimum c-value" ///				
					) ///
					modelwidth(30) note("Missing values imply one of the following: `=char(13)'`=char(10)'" ///
					"(1) That no original coefficient was estimated. `=char(13)'`=char(10)'" ///
					"(2) That no reversal is possible. `=char(13)'`=char(10)'" ///
					"(3) That the minimum reversing c-value falls outside the search range. `=char(13)'`=char(10)'")
				}
				else { // standard layout
					esttab matrix(`_to_display_gamma', fmt(3)), mtitles("") ///
					varlabels( ///
					1 "Untransformed gamma-shifted coeff.:" ///
					2 "Minimum cost:" ///				
					3 "Smallest c-value for reversal:"	///			
					) ///
					varwidth(35) note("Missing values imply one of the following: `=char(13)'`=char(10)'" ///
					"(1) That no original coefficient was estimated. `=char(13)'`=char(10)'" ///
					"(2) That no reversal is possible. `=char(13)'`=char(10)'" ///
					"(3) That the minimum reversing c-value falls outside the search range. `=char(13)'`=char(10)'")
				}
			}
		}	
	}

	
	
	*=====================================
	*16. Display to user if python is NOT specified
	*====================================
	
	if "`python'" != "python" {
		*-------------------------------------
		*16.1 Main display
		*-------------------------------------
		
		if "`fast'"=="" {
			noi {
				dis ""
				dis ""		
				dis "{bf:Results:}"
				if "`transpose'" != "" { // sets alternative layout suggested by Anthony
					esttab matrix(`_to_display', fmt(3) transpose), mtitles("") ///
					collabels( ///
					"Original coefficient" ///
					"Minimum c-value"  ///
					"Original p-value" ///
					"Maximum p-value" ///
					"Minimum p-value" ///
					"Min. c-value for insig." ///
					"Min. c-value for sig." ) ///
					modelwidth(25) note("Missing values imply one of the following: `=char(13)'`=char(10)'" ///
					"(1) That no original coefficient was estimated. `=char(13)'`=char(10)'" ///
					"(2) That no reversal is possible. `=char(13)'`=char(10)'" ///
					"(3) That the minimum reversing c-value falls outside the search range. `=char(13)'`=char(10)'" ///
					"(4) That coefficients cannot be made significant or insignificant within the search range.") 
				}
				else { // standard layout
					esttab matrix(`_to_display', fmt(3)), mtitles("") ///
					varlabels( ///
					1 "Original coefficient:" ///
					2 "Smallest c-value for reversal:"  ///
					3 "Original p-value:" ///
					4 "Maximum p-value:" ///
					5 "Minimum p-value:" ///
					6 "Smallest c-value for insig. coeff.:" ///
					7 "Smallest c-value for sig. coeff.:" ) ///
					varwidth(35) note("Missing values imply one of the following: `=char(13)'`=char(10)'" ///
					"(1) That no original coefficient was estimated. `=char(13)'`=char(10)'" ///
					"(2) That no reversal is possible. `=char(13)'`=char(10)'" ///
					"(3) That the minimum reversing c-value falls outside the search range. `=char(13)'`=char(10)'" ///
					"(4) That coefficients cannot be made significant or insignificant within the search range.") 
				}
			}
		}
		else {
			noi {
				dis ""
				dis ""		
				dis "{bf:Results:}"
				if "`transpose'" != "" { // sets alternative layout suggested by Anthony
					esttab matrix(`_to_display', fmt(3) transpose), mtitles("") ///
					collabels( ///
					"Original coefficient" ///
					"Minimum c-value"  ///
					"Original p-value" ///
					"Maximum p-value" ///
					"Minimum p-value") ///
					modelwidth(25) note("Missing values imply one of the following: `=char(13)'`=char(10)'" ///
					"(1) That no original coefficient was estimated. `=char(13)'`=char(10)'" ///
					"(2) That no reversal is possible. `=char(13)'`=char(10)'" ///
					"(3) That the minimum reversing c-value falls outside the search range. `=char(13)'`=char(10)'")
				}
				else { // standard layout
					esttab matrix(`_to_display', fmt(3)), mtitles("") ///
					varlabels( ///
					1 "Original coefficient:" ///
					2 "Smallest c-value for reversal:"  ///
					3 "Original p-value:" ///
					4 "Maximum p-value:" ///
					5 "Minimum p-value:") ///
					varwidth(35) note("Missing values imply one of the following: `=char(13)'`=char(10)'" ///
					"(1) That no original coefficient was estimated. `=char(13)'`=char(10)'" ///
					"(2) That no reversal is possible. `=char(13)'`=char(10)'" ///
					"(3) That the minimum reversing c-value falls outside the search range. `=char(13)'`=char(10)'" ///
					"(4) That coefficients cannot be made significant or insignificant within the search range.") 
				}
			}		
		}
		
		*-------------------------------------
		*16.2 Additional gamma display
		*-------------------------------------
		
		if `gamma' != 0 {
			noi {
				dis ""
				dis ""		
				dis "{bf:Further results after applying shift by gamma:}"
				if "`transpose'" != "" { // sets alternative layout suggested by Anthony
					esttab matrix(`_to_display_gamma', fmt(3) transpose), mtitles("") ///
					collabels( ///
					"Untrans. gamma-shifted coeff." ///
					"Minimum c-value" ///				
					) ///
					modelwidth(30) note("Missing values imply one of the following: `=char(13)'`=char(10)'" ///
					"(1) That no original coefficient was estimated. `=char(13)'`=char(10)'" ///
					"(2) That no reversal is possible. `=char(13)'`=char(10)'" ///
					"(3) That the minimum reversing c-value falls outside the search range. `=char(13)'`=char(10)'")
				}
				else { // standard layout
					esttab matrix(`_to_display_gamma', fmt(3)), mtitles("") ///
					varlabels( ///
					1 "Untransformed gamma-shifted coeff.:" ///
					2 "Smallest c-value for reversal:"	///			
					) ///
					varwidth(35) note("Missing values imply one of the following: `=char(13)'`=char(10)'" ///
					"(1) That no original coefficient was estimated. `=char(13)'`=char(10)'" ///
					"(2) That no reversal is possible. `=char(13)'`=char(10)'" ///
					"(3) That the minimum reversing c-value falls outside the search range. `=char(13)'`=char(10)'")
				}
			}
		}	
	}	

	*=====================================
	*17. Return results in r()
	*=====================================
	
	return matrix result `_to_display'
	if `gamma' != 0 return matrix resultgamma `_to_display_gamma'
	return matrix b `b_result' 
	return matrix r `r_result' 
	if "`dq'"=="" return matrix rgamma `r_gamma_result' 
	
	return matrix d `d_result' 
	return matrix hdp `hd_p_vals' 
	
	if "`fast'"=="" return matrix p `p_result' 
	if "`fast'"=="" return matrix y `y_result' 
	if "`fast'"=="" return matrix ygamma `y_gamma_result' 
	
	if "`python'"!="" return matrix cost 	`cost_mat'
	if "`python'"!="" return matrix labels 	`label_mat'
	
	return matrix bgamma `b_gamma_result' 
	if "`fast'"=="" return matrix pgamma `p_gamma_result' 
		
	*=====================================
	*18. Restore the original model.
	*=====================================

	qui estimates restore `prevmodel'
	
	}	// ends the qui condition
	
end
