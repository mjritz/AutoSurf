#!/usr/bin/env python
import sys
import numpy as np
import variables
import os
import subprocess
import random
from variables import *
#LAMMPS input script to automate simulations 
#Updated version should have less hard scripted information in script (i.e. pair coefficients, end=polar component, head = polar component, tail = alkyl group, surface =      surface)

def LAMMPS_long_equilibrate_monomer_surface(temp, z_distance,epsilon,sigma,r_distance, equil_time_long):
    
    #Generate a random number for the LAMMPS velocity input
    seed = random.randrange(1000000,9999999)

    h =""" 

variable T equal %6.2f

processors 2 2 *
units           real
boundary        p p p
atom_style      full
bond_style      harmonic
angle_style     harmonic

pair_style mie/cut  25.0
special_bonds lj 0.0 1.0 1.0

read_data system_conf.data

pair_coeff 1 1 0.68442  4.3635 15    6 25.0
pair_coeff 1 2 0.74619  4.0440 16.86 6 25.0
pair_coeff 1 3 0.496761 4.0547 10.75 6 25.0
pair_coeff 2 2 0.916247 3.7244 19    6 25.0
pair_coeff 2 3 0.977627 3.7352 11.94 6 25.0
pair_coeff 3 3 0.794739 3.7459 8     6 25.0
pair_coeff 3 4 0.95384  3.9036 11.94 6 25.00
pair_coeff 1 4 0.69948  4.2124 16.86 6 25.00
pair_coeff 2 4 0.77897  3.8929 19.00 6 25.00
pair_coeff 4 4 0.78871  4.0613 19.00 6 25.00

min_style       cg
neigh_modify    delay 10 every 1 check yes

group surface type 3
group nonpolar type 1
group terminal_polar type 2
group polar type 4

minimize 1.0e-6 1.0e-8 250 2508

fix 1 all nvt temp ${T} ${T} 100.0
velocity all create ${T} %i mom yes rot yes dist gaussian
velocity all zero linear
velocity all zero angular
fix zwall all wall/harmonic zlo -%4.2f %3.1f %3.1f %3.1f zhi %4.2f %3.1f %3.1f %3.1f pbc yes

thermo      %i
thermo_style custom step spcpu cpuremain pxx pyy pzz

timestep 2.0

dump  4 all dcd %i equil_initial.dcd
dump_modify 4 unwrap yes

fix 30 all balance 10000 1.0 shift z 100 1.0

fix 200 all momentum 100 linear 1 1 1

run %i
 
write_data system_conf.data
    """ % (temp, seed, float(z_distance),float(epsilon), float(sigma),float(r_distance), float(z_distance),float(epsilon),float(sigma),float(r_distance), int(0.05*equil_time_long), int(0.10*equil_time_long), equil_time_long)
    return h


def LAMMPS_short_equilibrate_monomer_surface(temp, z_distance,epsilon,sigma,r_distance, equil_dcd,equil_time_short):
    
    #Generate a random number for the LAMMPS velocity input
    seed = random.randrange(1000000,9999999)
    
    #LAMMPS input script to automate simulations 
    #Updated version should have less hard scripted information in script (i.e. pair coefficients, terminal=polar component, head = polar component, tail = alkyl group, surface = surface)
    
    h =""" 
variable T equal %6.2f

processors 2 2 *
units           real
boundary        p p p
atom_style      full
bond_style      harmonic
angle_style     harmonic

pair_style mie/cut  25.0
special_bonds lj 0.0 1.0 1.0

read_data system_conf.data

pair_coeff 1 1 0.68442  4.3635 15    6 25.0
pair_coeff 1 2 0.74619  4.0440 16.86 6 25.0
pair_coeff 1 3 0.496761 4.0547 10.75 6 25.0
pair_coeff 2 2 0.916247 3.7244 19    6 25.0
pair_coeff 2 3 0.977627 3.7352 11.94 6 25.0
pair_coeff 3 3 0.794739 3.7459 8     6 25.0
pair_coeff 3 4 0.95384  3.9036 11.94 6 25.00
pair_coeff 1 4 0.69948  4.2124 16.86 6 25.00
pair_coeff 2 4 0.77897  3.8929 19.00 6 25.00
pair_coeff 4 4 0.78871  4.0613 19.00 6 25.00

min_style       cg
neigh_modify    delay 10 every 1 check yes

group surface type 3
group nonpolar type 1 
group terminal_polar type 2
group polar type  4

#Next 5 lines "minimization" step
fix LIMIT all nve/limit 0.1
fix RESCALE all temp/rescale 10 ${T} ${T} 10.0 0.9

run 10000 

unfix LIMIT
unfix RESCALE

fix 1 all nvt temp ${T} ${T} 100.0
velocity all create ${T} %i mom yes rot yes dist gaussian
velocity all zero linear
velocity all zero angular
    """ % (temp, seed)
    
#Equilibrates system without an energy barrier blocking monomers from crossing the periodic boundary condition. 
    h+="""
thermo      %i
thermo_style custom step spcpu cpuremain temp etotal pe pxx pyy pzz

timestep 2.0

dump 4 all dcd %i %s
dump_modify 4 unwrap yes

fix 30 all balance 10000 1.0 shift z 100 1.0
fix 200 all momentum 100 linear 1 1 1

run %i 
    """ %(int(0.05*equil_time_short), int(0.10*equil_time_short), equil_dcd, equil_time_short)
    return h

def LAMMPS_density_profile(steps, RDF, prod_dcd):
    #Generates data to analyze the density distribution of the three components (3=surface, 2=polar component, 1=nonpolar component)
    h="""
undump 4
unfix 200
unfix 30

# reset timestep to zero so fix averaging works properly!  <3
reset_timestep 0

variable production equal %s
fix MOM all momentum 10 linear 0 0 1

variable   Nevery              equal 20
variable   Nfreq               equal ${production}-10*${Nevery}
variable   Nrepeat             equal floor(${Nfreq}/${Nevery})-1

compute cc1 terminal_polar chunk/atom bin/1d z lower 0.25
fix cc2 terminal_polar ave/chunk ${Nevery} ${Nrepeat} ${Nfreq} cc1 density/number file terminal_polar.dat

compute cc3 nonpolar chunk/atom bin/1d z lower 0.25
fix cc4 nonpolar ave/chunk ${Nevery} ${Nrepeat} ${Nfreq} cc3 density/number file nonpolar.dat

compute cc5 surface chunk/atom bin/1d z lower 0.25
fix cc6 surface ave/chunk ${Nevery} ${Nrepeat} ${Nfreq} cc5 density/number file surface.dat

compute NEIGHBOR terminal_polar coord/atom %3.1f 3
fix HW terminal_polar ave/histo ${Nevery} ${Nrepeat} ${Nfreq} -0.5 50.5 50 c_NEIGHBOR file histogram.dat mode vector

thermo          %i
thermo_style custom step spcpu cpuremain temp etotal pe pxx pyy pzz 

timestep        2.0

dump 10 all dcd %i %s
dump_modify 10 unwrap yes

#run balance command every 10000 steps in z direction
fix 20 all balance 1000 1.0 shift z 100 1.0

run ${production} 

variable h0 equal f_HW[1][3] # 0 bin of histogram
variable h1 equal f_HW[2][3] # 1 bin of histogram 
variable h2 equal f_HW[3][3] # 2 bin of histogram
variable h3 equal f_HW[4][3] # 3 bin of histogram
variable h4 equal f_HW[5][3] # 4 bin of histogram
variable h5 equal f_HW[6][3] # 5 bin of histogram

variable N equal (count(terminal_polar))/2 #divide by number of end groups per monomer and for two surfaces

print "$N ${h0} ${h1} ${h2} ${h3} ${h4} ${h5}" append hist_data.dat

write_data system_conf.data
    """ % (steps, RDF,(0.05*steps), (0.10*steps), prod_dcd)
    return h

def ST_input_script(z_distance, epsilon, sigma, r_distance, equil_dcd, equil_time_short, bead_per, prod_dcd, per_loop, steps): 
    string ="%d %g %20.15g %20.15g %20.15g" 
    
    h = """

thermo      %i
thermo_style custom step spcpu cpuremain etotal ke temp  pe press  pxx pyy pzz

timestep 2.0

dump  4 all dcd %i %s
dump_modify 4 unwrap yes

#run balance command every 10000 steps in z direction
fix 30 all balance 10000 z 100 1.0
fix 200 all momentum 100 linear 1 1 1

run %i

write_data minimized.data
    """ %(int(0.05*equil_time_short), int(0.10*equil_time_short), equil_dcd,equil_time_short)
    
    h+="""
unfix 200

reset_timestep 0
variable input_name string system
variable mol_size string %s
variable delta string 0.0005
variable xyz_file string xyz.txt
variable xyz_file_surface string xyz_surface.txt
variable xyz_p_file string xyz_p.txt
variable xyz_m_file string xyz_m.txt

neigh_modify    delay 10 every 1 check yes

group surface type 3
group monomer type 1 2 4

minimize 1.0e-6 1.0e-8 20 20

thermo          0
thermo_style custom step spcpu etotal ke temp pe press pxx pyy pzz cella cellb cellc

dump 10 all dcd %i %s
dump_modify 10 unwrap yes

#open dE.txt file and write header
print "#step pxx pyy pzz" file P.txt screen no
print "#step U Up Um U" file dE.txt screen no
print "#step volume_fraction" file vol_frac.txt screen no
run 0

#define variables for half box lengths
variable x equal cella/2.0
variable boxx2 equal $x
variable x equal cellb/2.0
variable boxy2 equal $x
variable x equal cellc/2.0
variable boxz2 equal $x
variable x equal cella
variable boxx equal $x
variable x equal cellb
variable boxy equal $x
variable x equal cellc
variable boxz equal $x
#define variables for plus and minus perturbation amounts
variable pd equal sqrt(1+${delta})
variable md equal sqrt(1-${delta})

#new
variable pvecx equal pxx
variable pvecy equal pyy
variable pvecz equal pzz
variable i equal step
#--

#number of TA evaluations to run
variable a loop %s
label loop

#number of steps to run between TA evaluations
run %s

#define current step
variable x equal step
variable frame equal $x

#define current potential energy without perturbation
#new
print "${i} ${pvecx} ${pvecy} ${pvecz}" append P.txt screen no
#--
variable x equal pe
variable U equal $x

#write out unscaled monomer coordinates to file name defined by variable xyz_file
dump 100 monomer custom 1000 ${xyz_file} id mass xu yu zu
dump 101 surface custom 1000 ${xyz_file_surface} id mass xu yu zu
dump_modify 100 sort id first yes format "%s"
dump_modify 101 sort id first yes format "%s"
run 0
undump 100
undump 101

#use python script to rescale monomer molecule coordinates
shell python rescale.py ${xyz_file} ${xyz_p_file} ${xyz_m_file} ${delta} ${mol_size} ${boxx} ${boxy} ${boxz}

#rescale surface coordinates
change_box surface x scale ${pd} z volume y scale ${pd} z volume remap
#read in rescaled  coordinates
read_dump ${xyz_p_file} 0 x y z box no replace yes scaled no wrapped no
run 0
variable x equal pe
variable Up equal $x

#scale surface coordinates to minus perturbation (this could be done in 1 step)
change_box surface x final -${boxx2} ${boxx2} y final -${boxy2} ${boxy2} z final -${boxz2} ${boxz2} remap
change_box surface x scale ${md} z volume y scale ${md} z volume remap
#read in rescaled  coordinates
read_dump ${xyz_m_file} 0 x y z box no replace yes scaled no wrapped no
run 0
variable x equal pe
variable Um equal $x

#reset box 
change_box surface x final -${boxx2} ${boxx2} y final -${boxy2} ${boxy2} z final -${boxz2} ${boxz2} remap
read_dump ${xyz_file} ${frame} x y z box no replace yes scaled no wrapped no

#calculate potential energy after
run 0
variable x equal pe
variable U2 equal $x 
print "${frame} ${U} ${Up} ${Um} ${U2}" append dE.txt screen no
variable frac file temp_vol_frac.txt
print "${frame} ${frac}" append vol_frac.txt screen no
variable frac delete
next a
jump SELF loop
    """ % (bead_per, int(0.10*steps), prod_dcd, per_loop, steps, string, string)  
    return h

RDF, xbox, ybox, zbox, surface_filename, surface_mol, pack_tol, bead_dist, boxx, boxy, boxz, bead_per, main_branch, equil_time_short, equil_time_long,monomer_x, monomer_y, monomer_zlow,          monomer_zhi,per_loop, steps,temp, z_distance, epsilon, sigma, r_distance, equil_dcd, prod_dcd, cutoff = variables.read_variables('input/input_file.txt')

Equil=LAMMPS_long_equilibrate_monomer_surface(temp, z_distance,epsilon,sigma,r_distance, equil_time_long)

STscript=LAMMPS_short_equilibrate_monomer_surface(temp, z_distance,epsilon,sigma,r_distance,equil_dcd,  equil_time_short)
STscript+=ST_input_script(z_distance, epsilon, sigma, r_distance, equil_dcd, equil_time_short, bead_per, prod_dcd, per_loop, steps)

Dscript=LAMMPS_short_equilibrate_monomer_surface(temp, z_distance,epsilon,sigma,r_distance,  equil_dcd,equil_time_short)
Dscript+=LAMMPS_density_profile(steps, RDF, prod_dcd)

#Read in the  structure and input file for moltemplate 
Density_lammps_file = open('system_density.in', 'wa')

#Read in the  structure and input file for moltemplate 
Initial_Equil_file = open('system_initial.in', 'wa')


#Read in the  structure and input file for moltemplate 
ST_lammps_file = open('system_ST.in', 'wa')

Initial_Equil_file.write(Equil)
Initial_Equil_file.close()
Density_lammps_file.write(Dscript)
Density_lammps_file.close()
ST_lammps_file.write(STscript)
ST_lammps_file.close()
