#Attempt to find interface using each frame from dcd file 

#!/usr/bin/env python 

import sys
import numpy as np
import ReadDCD
import variables

f2 = open('monomer_surface.xyz', 'wa') #Open in write and append format
outstream = f2

#Read in dcd file, define local surface, calculate distance from monomers
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
#print monomer_mol
for monomer_mol in range(surface_mol, total_molecules):
        monomers.append(xyz_final_frame[monomer_mol,:])
        monomer_mol+=1
        monomer_array =np.array(monomers, float)

#FOR DEBUGGING

#print xyz_final_frame
#print xyz_final_frame.shape
#print surface_array
#print surface_array.shape

#create array for surface molecules and monomer molecules 
#x = []
#y = []
#z = []
#x.append(xyz_final_frame[:,0])
#y.append(xyz_final_frame[:,1])
#z.append(xyz_final_frame[:,2])
#x_array = np.array(x, float)
#y_array = np.array(y, float)
#z_array = np.array(z, float)
#print x_array.shape


#concatenate the four [N,1] arrays together to form and [N,4] array 
#xyz_all_array = np.concatenate((x_array, y_array, z_array), axis=0)

#print xyz_all_array.shape
#print xyz_all_array

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

    #calculate distribution of surfaces as a function of z
    delta_z =0.50
    surface_freq,bins=np.histogram(xyz_surface[:,2],bins=np.linspace(-float(boxD[2]/2.0), float(boxD[2]/2.0)+1, ((float(boxD[2]/2.0)+1)-(-float(boxD[2]/2.0)))/delta_z))
    #print surface_freq

    #find min and max of surface, cutoff used arbitrarily to avoid pure   vapor
    #assumes that there is a value of more than 20 surfaces in an 80x80x1 A slab (from ND, fine assumption), 20/(80x80x1) = 0.003125
    #Using density, NA, g/mol, 1bead/molecule => 1bead/59.88Ang => should be more than 30 molecules in 60x60x1 slab, 30/(60x60x1)
    vapor_cut=.01388*boxD[0]*boxD[1]
    minl=np.argmax(surface_freq>vapor_cut)
    minr=np.size(surface_freq)-np.argmax(surface_freq[::-1]>vapor_cut)
    #print minl, minr

    #estimate bulk surface number density from middle 20 A of surface slab
    slabcenter=(minl+minr)/2.0
    surface_density=np.average(surface_freq[slabcenter-10:slabcenter+10])
    #surface density should be around 60 for 60x60x1 slab 
    #print surface_density

    #find interface
    interface_l=np.argmax(surface_freq>surface_density/2.0)
    interface_r=np.size(surface_freq)-np.argmax(surface_freq[::-1]>surface_density/2.0)-1
    #print interface_l, interface_r

    lower_interface = -int(boxD[2]/2.0)+(delta_z*interface_l)
    upper_interface = -int(boxD[2]/2.0)+(delta_z*interface_r)
    return lower_interface, upper_interface    

lower_interface, upper_interface = find_interface(boxx, boxy, boxz, surface_array) 

top=[]
bottom=[]
num_monomer_mol = monomer_array.shape[0]
for zvalue in monomer_array[:,2]:
    if zvalue>0:
        top.append(zvalue)
    if zvalue<0:
        bottom.append(zvalue)
total_top = len(top)/bead_per_monomer
total_bottom = len(bottom)/bead_per_monomer
print total_top

num_monomer_mol = monomer_array.shape[0]/bead_per_monomer
pdb_array = np.concatenate((surface_array,monomer_array), axis=0)
number_of_atoms = str(pdb_array.shape[0])
comment = "Created after Adding monomers"

#Writing out xyz file -----------------------------------------------------
s=number_of_atoms
s+= "\n"
s+=comment
s+="\n"
for (i, (x, y, z)) in enumerate(pdb_array):
    s+="%3i%12.3f%12.3f%12.3f\n" % (i+1, x, y, z)

outstream.write(s)
outstream.close()

