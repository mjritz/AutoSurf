#!/usr/bin/env python

# Adapted from script: press2surfTens.py by Mario Orsi (m.orsi at qmul.ac.uk, www.orsi.sems.qmul.ac.uk)
# Purpose: Calculates the surface tension using the standard formula:
#          0.5*Lz*(Pzz-0.5(Pxx+Pyy))
# Syntax: 
# Example: 
# Notes: Pxx, Pyy and Pzz are provided in units of atm (LAMMPS
#        convention for 'units real' style) and Lz in Angstrom
# Reference: Ismail et al, J Chem Phys 125, 014702 (2006) 

import sys, string
from math import sqrt
import numpy as np

P_data=np.loadtxt('P.txt')
print "Insert Lz"
Lz = str(raw_input('> '))

atm_in_Pa = float(101325) # note: 1 Pa = 1 N/m^2
A_in_m = float(1e-10) # Angstrom in meter
N_in_mN = float(1e3) # Newton in milliNewton

Lz = float(Lz) * A_in_m
print P_data[:,1]
Pxx = (P_data[:,1]) #* atm_in_Pa 
Pyy = (P_data[:,2]) #* atm_in_Pa 
Pzz = (P_data[:,3]) #* atm_in_Pa
Pxx_avg = np.mean(Pxx) * atm_in_Pa
Pyy_avg = np.mean(Pyy) * atm_in_Pa
Pzz_avg = np.mean(Pzz) * atm_in_Pa

surfTens = 0.5 * Lz * ( Pzz_avg - 0.5*(Pxx_avg+Pyy_avg) ) # [N/m]
surfTens = surfTens * N_in_mN # [mN/m]

print "Surface tension: % .6f mN/m" % surfTens

