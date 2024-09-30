#*******************************************************************************
#Reversing the reversal
#*******************************************************************************
#Python parts for reversal utility 
#*******************************************************************************

#=====================================
#1. Set-up
#=====================================

#-------------------------------------
#1.1 Import libraries
#-------------------------------------

import __main__
import os
os.environ["KMP_DUPLICATE_LIB_OK"]="TRUE"
import pandas as pd
import numpy as np
from sfi import Data
from sfi import Macro
from sfi import Matrix

from scipy.optimize import minimize
from scipy.optimize import LinearConstraint
from scipy.optimize import NonlinearConstraint
from scipy.optimize import BFGS

#============================================================
# 2. Define cost function
#============================================================

cost_type=int(Macro.getLocal('cost_type'))

#-------------------------------------
#1.1 Define cost function (using squared differences)
#-------------------------------------

if cost_type==1:
    def cost(l):
        
        #-------------------------------------
        #1.2 Basics
        #-------------------------------------
        K = len(l)
        minl = l[0]
        maxl = l[K-1]
        
        #-------------------------------------
        #1.3 dl
        #-------------------------------------
        dl = [0]*(K-1)
        for i in range(0,K-1):
            dl[i] = l[i+1] - l[i]
        maxdl = maxl - minl
        N = K - 1
        maxvar = (1/N - 1/N**2)*maxdl**2
        Edl = maxdl/N

        #-------------------------------------
        #1.4 Compute each component
        #-------------------------------------
        var_comp = [0]*(N)
        for i in range(0,N):
            var_comp[i] = (dl[i] - Edl)**2 
    
        #-------------------------------------
        #1.5 Get the variance, normalise, and output
        #-------------------------------------
        var = (1/N)*np.sum(var_comp)
        cost = var/maxvar
        return cost

#-------------------------------------
#1.2 Define cost function (using absolute value of differences)
#-------------------------------------

if cost_type == 2:
    def cost(l):
        
        #-------------------------------------
        #1.2 Basics
        #-------------------------------------
        K = len(l)
        minl = l[0]
        maxl = l[K-1]
        
        #-------------------------------------
        #1.3 dl
        #-------------------------------------
        dl = [0]*(K-1)
        for i in range(0,K-1):
            dl[i] = l[i+1] - l[i]
        maxdl = maxl - minl
        N = K - 1
        maxvar = (1/N - 1/N**2)*maxdl**2
        Edl = maxdl/N

        #-------------------------------------
        #1.4 Compute each component
        #-------------------------------------
        var_comp = [0]*(N)
        for i in range(0,N):
            var_comp[i] = (dl[i] - Edl)**2 
    
        #-------------------------------------
        #1.5 Get the variance, normalise, and output
        #-------------------------------------
        var = (1/N)*np.sum(var_comp)
        cost = (var/maxvar)**0.5
        return cost

#=====================================
#3. Definitions and imports
#=====================================

#-------------------------------------
#3.1 Import coefficients, scale_min, and scale_max
#-------------------------------------

#Coefficients
bdm = np.asarray(Matrix.get("_numerator"))[0:]
bdm = bdm.tolist()
bdm = [val[0] for val in bdm]

bdn = np.asarray(Matrix.get("_denominator"))[0:]
bdn = bdn.tolist()
bdn = [val[0] for val in bdn]

bdm = np.asarray(bdm)
bdn = np.asarray(bdn)

#Number of labels
nlabs = len(bdm)+1

#Original labels
l_original    = np.array(range(1,nlabs+1,1))

#scale_min and scale_max
scale_min = np.amin(l_original)
scale_max = np.amax(l_original)

#Get ratio at cost=0
b_numer = (l_original[0]-l_original[1])*bdm[0]
b_denom = (l_original[0]-l_original[1])*bdn[0]
for i in range(1,nlabs-1):
	j = i+1
	b_numer += (l_original[i]-l_original[j])*bdm[i]
	b_denom += (l_original[i]-l_original[j])*bdn[i]
	ratio0 = b_numer/b_denom	

#Get original coefficients
def originator(x):
	original = (l_original[0]-l_original[1])*x[0]
	for i in range(1,nlabs-1):
		j = i+1
		original += (l_original[i]-l_original[j])*x[i]
	return original
bdm_original = originator(bdm)
bdn_original = originator(bdn)

#Get max and min ratio 
ratios = bdm/bdn
max_ratio = np.amax(ratios)
min_ratio = np.amin(ratios)

#Get precision 
precision = int(Macro.getLocal('precision'))

#-------------------------------------
#3.2 Set constraint that labels are weakly increasing
#-------------------------------------
	
#Define monotonicity constraints
#This will create a matrix with rows and columns equal to the number K of labels, with the first K-1 rows having elements 1 at k and -1 at k+1. The last row has all zeros. I think the last row is an error caused above, but then cancels out. But the thing works correctly for now...
monotone_array1 = []
x = [0]*nlabs
for i in range(0,nlabs-1):
	j = i+1
	y = x.copy()
	y[i] = 1
	y[j] = -1	
	monotone_array1.append(y) 
monotone_array1.append(x) 
monotone_array2 = [-np.inf]*nlabs
monotone_array3 = [0]*nlabs

monotonicity_constraint = LinearConstraint(monotone_array1, monotone_array2, monotone_array3)

#-------------------------------------
#3.3 Define constraint that new labels need to lead to a given coefficient ratio (set in the wrapper below)
#-------------------------------------

def ratio_cons_f(x):
	value_numer = (x[0]-x[1])*bdm[0]
	value_denom = (x[0]-x[1])*bdn[0]
	for i in range(1,nlabs-1):
		j = i+1
		value_numer += (x[i]-x[j])*bdm[i]
		value_denom += (x[i]-x[j])*bdn[i]
	value = value_numer/value_denom	
	return value	
	
#-------------------------------------
#3.4 Set constraint that labels need to be between scale_min and scale_max. 
#-------------------------------------

tmp1 = [0]*nlabs
tmp1[0]=1
tmp2 = [0]*nlabs
tmp2[nlabs-1]=1

boundary_constraint = LinearConstraint([tmp1,tmp2], [scale_min,scale_max], [scale_min,scale_max])

#-------------------------------------
#3.5 Make a wrapper function for minimisation so that we only have one variable argument.  
#-------------------------------------

def minimize_wrapper_min(x):
	ratio_constraint_nonlinear = NonlinearConstraint(ratio_cons_f, -np.inf, x, jac='2-point', hess=BFGS()) 
	result = minimize(cost, l_transformed, constraints=[monotonicity_constraint, ratio_constraint_nonlinear, boundary_constraint], tol=1e-8, options = {'maxiter': 10000, 'disp': False})
	return result

#-------------------------------------
#3.6 Make a wrapper function for maximisation so that we only have one variable argument.  
#-------------------------------------

def minimize_wrapper_max(x):
	ratio_constraint_nonlinear = NonlinearConstraint(ratio_cons_f, x, np.inf, jac='2-point', hess=BFGS()) 
	result = minimize(cost, l_transformed, constraints=[monotonicity_constraint, ratio_constraint_nonlinear, boundary_constraint], tol=1e-8, options = {'maxiter': 10000, 'disp': False})
	return result

#=====================================
#4. Do the computations 
#=====================================
	
#-------------------------------------
#4.1 Minimisation 
#-------------------------------------

#Get a list of ratios over which to minimise
ratio_list = list(np.linspace(min_ratio, ratio0, num=precision, endpoint=True))
cost_list = []

# Iterate over those ratios and minimize the cost
counter = 0
for i in ratio_list:
	if counter==0:
		# Set l_transformed
		l_transformed = np.random.uniform(low=l_original[0], high=l_original[nlabs-1], size=nlabs) 
		l_transformed = np.sort(l_transformed)
		l_transformed[0] = l_original[0]
		l_transformed[nlabs-1] = l_original[nlabs-1]
	elif counter>0:
		l_transformed = resulti.x
	
	resulti = minimize_wrapper_min(i)
	cost_list.append(resulti.fun) 
	
	counter += 1

result_min = np.dstack((cost_list, ratio_list))[0, :]
	
#Sort
ind=np.argsort(result_min[:,0])
result_min=result_min[ind]

#-------------------------------------
#4.2 Maximisation 
#-------------------------------------

#Get a list of ratios over which to minimise
ratio_list = list(np.linspace(max_ratio, ratio0, num=precision, endpoint=True))
cost_list = []

# Iterate over those ratios and minimize the cost
counter = 0
for i in ratio_list:
	if counter==0:
		# Set l_transformed
		l_transformed = np.random.uniform(low=l_original[0], high=l_original[nlabs-1], size=nlabs) 
		l_transformed = np.sort(l_transformed)
		l_transformed[0] = l_original[0]
		l_transformed[nlabs-1] = l_original[nlabs-1]
	elif counter>0:
		l_transformed = resulti.x
	
	resulti = minimize_wrapper_max(i)
	cost_list.append(resulti.fun) 
	
	counter += 1

result_max = np.dstack((cost_list, ratio_list))[0, :]
	
#Sort
ind=np.argsort(result_max[:,0])
result_max=result_max[ind]

#=====================================
#5. Send to Stata 
#=====================================

Matrix.store("result_min", result_min)
Matrix.store("result_max", result_max)

Macro.setLocal("ratio_0", str(ratio0))
Macro.setLocal("min_ratio", str(min_ratio))
Macro.setLocal("max_ratio", str(max_ratio))
