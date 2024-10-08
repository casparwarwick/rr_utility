#*******************************************************************************
#Reversing the reversal
#*******************************************************************************
#Python parts for reversal utility 
#*******************************************************************************

#-------------------------------------
#1. Set-up
#-------------------------------------

import __main__
import os
os.environ["KMP_DUPLICATE_LIB_OK"]="TRUE"
import warnings
import pandas as pd
import numpy as np
from sfi import Data
from sfi import Macro
from sfi import Matrix

from scipy.optimize import minimize
from scipy.optimize import LinearConstraint
from scipy.integrate import quad

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
	

#-------------------------------------
#1.6 Import coefficients, number of labels, gamma, scale_min and scale_max, and original labels from Stata
#-------------------------------------

#Coefficients
bd = Data.get(var=["bd"])
sign_var = Data.get(var=["sign"])
sign = sign_var[0]

#Number of labels
nlabs = len(bd)  

#Gamma
gam = - float(Macro.getLocal('gamma'))
isgamma = float(Macro.getLocal('is_this_gamma'))
if isgamma==0:
	gam = 0

#scale_min and scale_max
scale_min = float(Macro.getLocal('scale_min'))
scale_max = float(Macro.getLocal('scale_max'))

#Original labels
labels_depvar_rescaled = np.asarray(Matrix.get("_labels_depvar_rescaled"))[0:]
labels_depvar_rescaled = labels_depvar_rescaled.tolist()
l_original    = [val[0] for val in labels_depvar_rescaled]
l_transformed = [val[0] for val in labels_depvar_rescaled] # for an initial value

#-------------------------------------
#1.7 Set constraint that labels are weakly increasing
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
#1.8 Set constraint that new labels need to lead to a reversal
#-------------------------------------

reversal_array = [0]*nlabs
for i in range(0,nlabs):
	h = i - 1
	if i==0:
		reversal_array[i] = bd[i]	
	if i>0 & i<(nlabs-1):
		reversal_array[i] = bd[i]-bd[h]
	if i==(nlabs-1):
		reversal_array[i] = -bd[h]
		
if sign > 0:
	reversal_constraint = LinearConstraint(reversal_array, [-np.inf], [sign*gam])
else:
	reversal_constraint = LinearConstraint(reversal_array, [sign*gam], [np.inf])
	
#-------------------------------------
#1.9 Set constraint that labels need to be between 1 and the width of the scale, given by "width". 
#-------------------------------------

tmp1 = [0]*nlabs
tmp1[0]=1
tmp2 = [0]*nlabs
tmp2[nlabs-1]=1

boundary_constraint = LinearConstraint([tmp1,tmp2], [scale_min,scale_max], [scale_min,scale_max])

#-------------------------------------
#1.10 Minimize cost function subject to the constraints
#-------------------------------------

#Minimize cost function and save result
result = minimize(cost, l_transformed, constraints=[monotonicity_constraint, reversal_constraint, boundary_constraint])

#Save cost
cost = [result.fun]*nlabs # just puts things into the right format 			

#-------------------------------------
#1.11 Output result to Stata
#-------------------------------------

Data.addVarDouble("python_labels")
Data.addVarDouble("python_cost")

Data.store("python_labels", None, result.x, None)
Data.store("python_cost", None, cost, None)