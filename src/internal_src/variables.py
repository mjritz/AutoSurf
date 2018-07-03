#!/usr/bin/env python 

import struct
import numpy as np

def read_variables(filename):
    var_file=open(filename,"rb")
    for vars in var_file.readlines():
        if '#' not in vars:
            variable = vars.split(' ')
            if variable[0] == 'rdf':
                RDF = float(variable[1])
            if variable[0] == 'surface_boxsize':
                xbox = float(variable[1])
                ybox = float(variable[2])
                zbox = float(variable[3])
            if variable[0] == 'surface_filename':
                surface_filename = variable[1]
            if variable[0] == 'surface_molecules':
                surface_mol = int(variable[1])
            if variable[0] == 'packmol_tolerance':
                pack_tol = float(variable[1])
            if variable[0] == 'bead_distance':
                bead_dist = variable[1]
            if variable[0] == 'overall_system_size':
                boxx = float(variable[1])
                boxy = float(variable[2])
                boxz = float(variable[3])
            if variable[0] == 'units_per_molecule':
                bead_per = variable[1]
            if variable[0] == "units_in_linear_branch":
                main_branch = variable[1]
            if variable[0] == 'short_equilibrate':
                equil_time_short = variable[1]
            if variable[0] == 'long_equilibrate':
                equil_time_long = variable[1]
            if variable[0] == 'monomer_space':
                monomer_x = float(variable[1])
                monomer_y = float(variable[2])
                monomer_zlow = float(variable[3])
                monomer_zhi = float(variable[4])
            if variable[0] == 'number_ST_loops':
                per_loop = variable[1]
            if variable[0] == 'step_per_loop':
                steps = variable[1]
            if variable[0] == 'temperature':
                temp = float(variable[1])
            if variable[0] == 'wall':
                z_distance = variable[1]
                epsilon = variable[2]
                sigma = variable[3]
                r_distance = variable[4]
            if variable[0] == 'DCD_filename':
                equil_dcd = variable[1]
                prod_dcd = variable[2]
            if variable[0] == 'cutoff_distance':
                cutoff = variable[1]
    return RDF, xbox, ybox, zbox, surface_filename, surface_mol, pack_tol, bead_dist, boxx, boxy, boxz, bead_per, main_branch, equil_time_short, equil_time_long, monomer_x, monomer_y, monomer_zlow, monomer_zhi, per_loop, steps, temp, z_distance, epsilon, sigma, r_distance, equil_dcd, prod_dcd, cutoff
