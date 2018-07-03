#!/usr/bin/env python

#This script reads in a surfactant structure file and outputs a data file for
#LAMMPS: 'python execute.py surfactant_structure.py' 
#In the same directory must have surfactant structure, surface.lt, mie_ff.lt, 
#'surface_name'.xyz, pdb_create.py, xyz_coord_bonds.py 


#Assumptions: no dendrimerization, no extra lines at the end of the structure file, T (main tail), B (side branches), H (head groups)
#Name of surface  file and tolerance are written in pdb_create.py near the bottom
#Initial distance of beads written in xyz_coord_bonds.py
#Size of system is written near the bottom of pdb_create.py

import write_lammps_input
import xyz_coord_bonds
import os
from pdb_create import i
import sys
from xyz_coord_bonds import head_groups
from xyz_coord_bonds import main_tail
from xyz_coord_bonds import descript

version = '1.0'
description = descript 
length = main_tail
head = head_groups
home_directory = '../'

os.system('cp input/*.lt .' )

if description != 'linear':
    print "Description of Branched Surface Agent"
    further_description = str(raw_input('>'))
    
print "Subdirectory?"
lammps_input_files = str(raw_input('> '))

os.system('src/packmol/packmol < monomer_surface.inp')

os.system('src/moltemp_src.2013-3-03/moltemplate.sh  -nocheck -xyz monomer_surface.xyz  system_conf.lt')
os.system('rm src/*.pyc')

if description == 'branched':
    os.system('mkdir %s/%s' %(home_directory,description))
    os.system('mkdir %s/%s/%i%i' %(home_directory, description, head, length))
    os.system('mkdir %s/%s/%i%i/%s' %(home_directory,description, head, length,further_description))
    os.system('mkdir %s/%s/%i%i/%s/%i' %(home_directory,description, head, length,further_description, i))
    os.system('mkdir %s/%s/%i%i/%s/%i/%s' %(home_directory,description, head, length,further_description, i, lammps_input_files))
    os.system('cp -rf src %s/%s/%i%i/%s/%i/%s/.' %(home_directory,description,head, length, further_description, i,           lammps_input_files))
    os.system('cp -rf src/internal_src %s/%s/%i%i/%s/%i/%s/src/.' %(home_directory,description,head, length, further_description, i,           lammps_input_files))
    os.system('cp -rf input %s/%s/%i%i/%s/%i/%s/.' %(home_directory,description,head, length, further_description, i, lammps_input_files))
    os.system('cp -rf os_run_scripts %s/%s/%i%i/%s/%i/%s/.' %(home_directory,description,head, length, further_description, i, lammps_input_files))
    os.system('cp -rf analysis %s/%s/%i%i/%s/%i/%s/.' %(home_directory,description,head, length, further_description, i, lammps_input_files))
    os.system('rm -r system_conf.in*')
    os.system('cp * %s/%s/%i%i/%s/%i/%s/.' %(home_directory,description,head, length, further_description, i,lammps_input_files))
    os.system('rm -rf output_ttree')
    os.system('rm *')
    os.system('clear')
    os.chdir('%s/%s/%i%i/%s/%i/%s/' % (home_directory, description,head, length, further_description, i,lammps_input_files))
    os.system('bsub < os_run_scripts/run_lammps_add')

#if description == 'linear': 
