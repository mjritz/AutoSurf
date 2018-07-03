#RDF analysis outputs a plot that has no slope up until CSA and then has a constant downward slope - trying to find # surfactant at which the slope change occurs. 

#!/usr/bin/env python 

import sys
import numpy as np
import os
import scipy
from scipy import stats
from numpy import ndarray
import matplotlib.pyplot as plt

#Read in information from output file from RDF analysis of surfactant and assign data to numpy arrays:
histogram_instream = np.loadtxt("hist_data.dat", skiprows=1)

num_surf_descending=histogram_instream[:,0]
hist1_descending=histogram_instream[:,1]
hist2_descending=histogram_instream[:,2]

#DEBUG

#flip vector so that the order is increasing
num_surf_ascending= np.fliplr([num_surf_descending])[0]
hist1_ascending= np.fliplr([hist1_descending])[0]
hist2_ascending= np.fliplr([hist2_descending])[0]
num_inputs = np.size(num_surf_ascending)

#DEBUG
#print hist2_ascending

#Error Analysis Linear--------------------
hist_difference_lin =[]
hist_difference_lin = hist2_ascending-hist1_ascending
error_array_lin=np.zeros([num_inputs,2])
error_difference =[]

#Error Analysis Slope--------------------
hist_difference_slope =[]
hist_difference_slope = hist2_descending-hist1_descending
hist_addition=hist2_descending+hist1_descending
linear_regression=np.zeros([num_inputs,2])
error_array_slope=np.zeros([num_inputs-1,4])
slope_intercept_array = np.zeros([num_inputs-1,3])
slope = 0
intercept = 0

#DEBUG
#print hist_difference_lin

#print num_surf_descending
#print hist_difference_slope

def solve_for_y(poly_coeffs, y):
	pc = poly_coeffs.copy()
	pc[-1] -=y
	return np.roots(pc)

z = np.polyfit(num_surf_descending, hist_addition, 3)
print num_surf_descending, hist_addition
print '0.10'
#print np.roots(z)
print solve_for_y(z, 0.10)

#cycle through possible options for # of objects on each side of intersect:
for n in range(1,num_inputs+1):
    difference_slope = hist_difference_slope[0:n]
    difference_lin = hist_difference_lin[0:n]
    average_lin =np.average(difference_lin)
    linear_regression[n-1]=(float(num_surf_descending[n-1]),difference_slope[n-1])
    if n==1 or n==2:
        pass
    else:
        xbar = ( ndarray.transpose(linear_regression[0:n,[0]]))
        ybar = ndarray.transpose(linear_regression[0:n,[1]])
        xbar_x_ybar = ndarray.transpose(linear_regression[0:n,[0]])*ndarray.transpose(linear_regression[0:n,[1]])
        xbar_squared = ndarray.transpose(linear_regression[0:n,[0]])**2
        n_inverse = (n**-1)
        denominator = (np.sum(xbar_squared)-(n_inverse*(np.sum(xbar)**2)))
        numerator = (np.sum(xbar_x_ybar)-(n_inverse*np.sum(xbar)*np.sum(ybar))) 
        beta = numerator/denominator
        ymean= np.average(ybar)
        xmean= np.average(xbar)
        alpha = (ymean - beta*xmean)
        intercept = alpha
        slope = beta
        slope_intercept_array[n-2] = (n,slope,intercept)
    error_slope=((slope*ndarray.transpose((linear_regression[0:n,[0]])))+(intercept-difference_slope))**2
    error_difference= difference_lin-average_lin
    error_lin = (error_difference**2)
    error_array_lin[n-1] = (n+1,np.sum(error_lin))
    error_array_slope[n-2]=(n,np.sum(error_slope),slope,intercept)
    n+=1

#DEBUG
#print error_array_lin

#-----------------------------------------

#DEBUG
#print error_array_slope

#-----------------------------------------

n=0
combined_error=np.zeros([num_inputs-1,4])

#Find the minimum error associated with the sloped line:
slope_attributes = len(error_array_slope)
for n in range(1,slope_attributes+1):
    combined_error[n-1]=(n, error_array_slope[num_inputs-n-1,1]+error_array_lin[n-1,1], error_array_slope[num_inputs-n-1,2],error_array_slope[num_inputs-n-1,3])
    n+=1

n=0
for n in range(1,slope_attributes+1):
    if combined_error[n-1,1] == np.amin(combined_error[:,1]):
        optimal_attributes = combined_error[n-1,0]


CSA="Error in code"
n=0
for n in range(1,slope_attributes+1):
    if combined_error[n-1,0]==optimal_attributes:
        CSA = abs(combined_error[n-1,3]/combined_error[n-1,2])

print CSA

#DEBUG 
#print combined_error
print optimal_attributes

