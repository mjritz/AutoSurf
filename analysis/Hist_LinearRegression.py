import matplotlib
import numpy as np
import sys
from scipy import stats
import os


#Read in histogram data
hist_data = np.loadtxt('hist_data.dat',skiprows=0)
x = (hist_data[:,0])
y = ((hist_data[:,1]+hist_data[:,2])*100)
n = 3

# Polynomial Regression
def polyfit(x, y, degree, percent):
    results = {}
    coeffs = np.polyfit(x, y, degree)
    # Polynomial Coefficients
    # r-squared
    p = np.poly1d(coeffs)
    # fit values, and mean
    yhat = p(x)                      # or [p(z) for z in x]
    ybar = np.sum(y)/len(y)          # or sum(y)/len(y)
    ssreg = np.sum((yhat-ybar)**2)   # or sum([ (yihat - ybar)**2 for yihat in yhat])
    sstot = np.sum((y - ybar)**2)    # or sum([ (yi - ybar)**2 for yi in y])
    error = ssreg / sstot
    results['Rsquared for 2D polynomial fit'] = ssreg / sstot
    percent_value = []
    for x in range(30,200, 1):
        if percent<((coeffs[0]*(x**3))+(coeffs[1]*(x**2))+(coeffs[2]*x)+(coeffs[3]))<(percent+0.5):
            percent_value.append(x)
    return coeffs, results,error, percent_value

coeffs, results, error, percent_value = polyfit(x,y,n,10)
print 'number of surfactants that have 0/1 water unit in the first hydration shell for 10-10.5% of the simulation.'
print percent_value
print 'error associated with polynomial regression curve - %i degree' %(n) 
print error
print 'maximum percentage calculated'
print np.amax(y)
print error

if np.amax(y) > 11 and error>0.9:
    os.sys('python src/execute_ST.py input/structure_file %i' %(percent_value[0]))
