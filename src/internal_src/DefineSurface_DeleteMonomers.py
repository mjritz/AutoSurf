#Define the location of the interface using each frame from dcd file 
#Create array with location of each monomer from the surface and delete monomers that are X angstroms from the defined interface. 

#!/usr/bin/env python 

import sys
import numpy as np
import ReadDCD
import variables

f2 = open('trimmed_system.xyz', 'wa') #Open in write and append format
outstream = f2

#Read in dcd file, define local surface, calculate distance from monomer 
RDF, xbox, ybox, zbox, surface_filename, surface_mol, pack_tol, bead_dist, boxx, boxy, boxz,     bead_per, main_branch, equil_time_short,equil_time_long, monomer_x, monomer_y, monomer_zlow,          monomer_zhi,per_loop, steps,temp, z_distance, epsilon, sigma, r_distance, equil_dcd,prod_dcd, cutoff = variables.read_variables('input/input_file.txt')

nAt,BoxD_frames,BoxAngle_frames,xyz_frames = ReadDCD.read_DCD(sys.argv[1])

#define total number of molecules, number of surface molecules, and number of monomers
total_molecules = nAt
surface_molecules = int(surface_mol)
bead_per_monomer = int(bead_per)
bead_in_main_branch = int(main_branch)
monomer_molecules = (total_molecules-surface_molecules)/bead_per_monomer
xyz_final_frame = xyz_frames[-1]

surface = []
surface_molecule=0
for surface_molecule in range(surface_mol):
    surface.append(xyz_final_frame[surface_molecule,:])
    surface_molecule+=1
surface_array =np.array(surface, float)

monomers = []
monomer_mol= surface_mol+1
for monomer_mol in range(surface_mol, total_molecules):
        monomers.append(xyz_final_frame[monomer_mol,:])
        monomer_mol+=1
        monomer_array =np.array(monomers, float)

g = 'original xyz file\n' 
for (i,(x,y,z)) in enumerate(surface_array):
    g += "%5i %8.3f%8.3f%8.3f\n" % (i+1, x, y, z)
for (i,(x,y,z)) in enumerate(monomer_array):
    g += "%5i %8.3f%8.3f%8.3f\n" % (i+1, x, y, z)
g+= "End\n"

#out_stream.write(g)
#out_stream.close()

def find_interface(boxx, boxy, boxz, surface_array):

    def wrapcoord(xyzcom,boxD,box_center=[0,0,0]):
        #shift box to center
        xyzcom-=box_center
        xyzcom-=boxD*np.around(xyzcom/boxD)
        return xyzcom
    
    D=3
    boxD=np.zeros(D)
    boxD[0]=float(boxx)*2
    boxD[1]=float(boxy)*2
    boxD[2]=float(boxz)

    #wrap surface coordinates to box
    xyz_surface=wrapcoord(surface_array,boxD,box_center=[0,0,0])

    #calculate distribution of surface as a function of z
    delta_z =0.50
    surface_freq,bins=np.histogram(xyz_surface[:,2],bins=np.linspace(-float(boxD[2]/2.0), float(boxD[2]/2.0)+1, ((float(boxD[2]/2.0)+1)-(-float(boxD[2]/2.0)))/delta_z))

    #find min and max of surface, cutoff used arbitrarily to avoid pure   vapor
    #assumes that there is a value of more than 20 molecules in an 80x80x1 A slab (from ND, fine assumption), 20/(80x80x1) = 0.003125
    #Using density, NA, g/mol, 1bead/molecule => 1bead/59.88Ang => should be more than 30 molecules in 60x60x1 slab, 30/(60x60x1)
    vapor_cut=.01388*boxD[0]*boxD[1]
    minl=np.argmax(surface_freq>vapor_cut)
    minr=np.size(surface_freq)-np.argmax(surface_freq[::-1]>vapor_cut)

    #estimate bulk surface number density from middle 20 A of surface slab
    slabcenter=(minl+minr)/2.0
    surface_density=np.average(surface_freq[slabcenter-10:slabcenter+10])
    #surface density should be around 60 for 60x60x1 slab 

    #find interface
    interface_l=np.argmax(surface_freq>surface_density/2.0)
    interface_r=np.size(surface_freq)-np.argmax(surface_freq[::-1]>surface_density/2.0)-1
    #print interface_l, interface_r

    lower_interface = -int(boxD[2]/2.0)+(delta_z*interface_l)
    upper_interface = -int(boxD[2]/2.0)+(delta_z*interface_r)
    return lower_interface, upper_interface    

lower_interface, upper_interface = find_interface(boxx, boxy, boxz, surface_array) 

accepted_monomers = []
mol_array = monomers_array.reshape(-1, bead_per_monomer, 3)

#search through the first beads to find  z values inside  of the range and copy
#molecule to a new list 
acceptable_range = float(cutoff) 
maximum = upper_interface-5
minimum = lower_interface+5
new_mol_number_top = (monomer_molecules/2)-10
new_mol_number_bottom = (monomer_molecules)-10
imol = 1

for zval in  mol_array[:,0,2]:
    if zval>maximum and  0<imol<=new_mol_number_top:
        accepted_monomer.append(mol_array[imol,:,:])
    if zval<maximum and new_mol_number_top+10<=imol<=new_mol_number_bottom:
        accepted_monomer.append(mol_array[imol,:,:])
    else: 
        pass
    imol+=1

top=[]
bottom=[]
accepted_monomer_array = np.array(accepted_monomer, float)
num_monomer_mol = accepted_monomer_array.shape[0]
for zvalue in accepted_monomer_array[:,0,2]:
    if zvalue>0:
        top.append(zvalue)
    if zvalue<0:
        bottom.append(zvalue)
total_top = len(top)
total_bottom = len(bottom)
good_monomer_array = accepted_monomer_array.reshape(-1, 3)
deleted = monomer_molecules - good_monomer_array.shape[0]/bead_per_monomer
pdb_array = np.concatenate((surface_array,good_monomer_array), axis=0)
number_of_atoms = str(pdb_array.shape[0])
comment = "Created after trimming monomers"

#Writing out xyz file -----------------------------------------------------
s=number_of_atoms
s+= "\n"
s+=comment
s+="\n"
for (i, (x, y, z)) in enumerate(pdb_array):
    s+="%3i%12.3f%12.3f%12.3f\n" % (i+1, x, y, z)

outstream.write(s)
outstream.close()

