#!/usr/bin/env python

import  DefineSurface 
import  variables
import sys

#print "Input file for moltemplate"
in_lt_file = open(sys.argv[2], 'w')
add = int(sys.argv[3])

RDF, xbox, ybox, zbox, surface_filename, surface_mol, pack_tol, bead_dist, boxx, boxy, boxz,     bead_per, main_branch, equil_time_short,equil_time_long, monomer_x, monomer_y, monomer_zlow, monomer_zhi,per_loop, steps,temp, z_distance, epsilon, sigma, r_distance, equil_dcd,prod_dcd, cutoff = variables.read_variables('input/input_file.txt')

g = 'import "surface.lt" \n'
g+= 'import "monomer_struct.lt" \n'

g+='''
write_once("Data Boundary") {
  -%4.6f   %12.6f xlo xhi
    -%4.6f   %12.6f ylo yhi
    -%4.5f   %12.5f zlo zhi
} \n\n''' % (boxx, boxx, boxy, boxy, boxz, boxz)

g+= "surface = new Surface [%i] \n\n" % (surface_mol)

g+= "monomer = new Monomer [%i]" % (int(DefineSurface.num_monomer_mol)+add)

in_lt_file.write(g)
in_lt_file.close()

packmol_file = open('monomer_additional.inp', 'wa')

#writing out packmol input file 
h = "tolerance %f\n" % (pack_tol)                                                                                           
h+="filetype xyz\n"                                                                                             
h+="output monomer_additional.xyz\n\n"""                                                                                                         

h += "structure monomer_surface.xyz" 
h += "  number 1\n"   

h += "  fixed 0. 0. 0. 0. 0. 0.\n"                                                                              
h += "end structure\n\n"                                                                                                          

h += "structure monomer.xyz\n"                                                                                       
h += "  number %i\n" % (add/2)                                                                                             
h += "  inside box -%4.1f -%4.1f %4.1f %4.1f %4.1f %4.1f \n" % (monomer_x, monomer_y, monomer_zlow, monomer_x, monomer_y, monomer_zhi)                                                             
h += "  end atoms\n"  
h += "end structure\n\n"                                                                                            
                                                                                                                   
h += "structure monomer.xyz\n"                                                                                       
h += "  number %i\n" %(add/2)                                                                                             
h += "  inside box -%4.1f -%4.1f -%4.1f %4.1f %4.1f -%4.1f \n" % (monomer_x, monomer_y, monomer_zhi, monomer_x, monomer_y, monomer_zlow)                                                          
h += "  end atoms\n"                                                                                              
h += "end structure\n\n"                                                                                            

packmol_file.write(h)
packmol_file.close()
