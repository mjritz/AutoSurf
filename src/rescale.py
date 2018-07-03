import numpy as np
import sys
import time
import bisect
import os

#this script assumes the following:
# 1 type of molecule
# molecules are in a continuous block
# assumes box center at 0,0,0
# assumes reading unwrapped coordinates
# Hard coded in: sigma of head group (line 167) and vapor pressure of surface (line 108) 

def mol_trans(xyz,delta,mol_com):
# translate com of molecules
	xyz_s=np.empty_like(xyz)
	scale=np.array([np.sqrt(1.0+delta),np.sqrt(1.0+delta),1.0/(1.0+delta)])
	
	transcom=(mol_com*scale)-mol_com

	for i in range(np.size(xyz,axis=0)):
		xyz_s[i]=xyz[i]+transcom[i]
	return np.array(xyz_s.reshape((-1,3)))

def write_lammps_dump(xyz,lammps_dump,header,dump_file_name):
	s=header
	for i,(pos) in enumerate(xyz):
		s+= ('%d %20.15g %20.15g %20.15g\n') % (lammps_dump[i,0],pos[0],pos[1],pos[2])
	f=open(dump_file_name,'w')
	f.write(s)
	f.close()

def read_lammps_dump(dump_file_name):
	f=open(dump_file_name,'r')
	dum=f.readlines()
	lfile=[]
	for i in dum:
		lfile.append(i.split())
	lammps_dump=np.array(lfile[9:],float)
	xyz=lammps_dump[:,2:5]
	mass=lammps_dump[:,1]
	return xyz,mass,lammps_dump

#wrap coordinates
def wrapcoord(xyzcom,boxD,box_center=[0,0,0]):
	#shift box to center
	xyzcom-=box_center
	xyzcom-=boxD*np.around(xyzcom/boxD)
	return xyzcom

def main():
	D=3
	surface_cut=40.0 #Changed this from 15!
	boxD=np.zeros(D)
	
	#set variables from command line inputs
	dump_file=sys.argv[1]
	xyz_p_file=sys.argv[2]
	xyz_m_file=sys.argv[3]
	delta=float(sys.argv[4])
	mol_size=int(sys.argv[5])
	boxD[0]=float(sys.argv[6])
	boxD[1]=float(sys.argv[7])
	boxD[2]=float(sys.argv[8])
	
	#read lammps dump file of monomer molecules
	xyz,mass,lammps_dump=read_lammps_dump(dump_file)
	
	#make lammps header for output. Everything except number of atoms and ITEM ATOMS has arbitrary values
	lammps_dump_header=("ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n"+str(np.size(lammps_dump,axis=0))+"\nITEM: BOX BOUNDS pp pp pp\n"
	"-40 40\n-40 40\n-150 150\nITEM: ATOMS id xu yu zu\n")

	#compute center of mass of each molecule
	mxyz=np.zeros_like(xyz)
	mol_mass=np.sum(mass[0:mol_size])
	for k in range(3):
		mxyz[:,k]=xyz[:,k]*mass
	mxyz=mxyz.reshape((-1,mol_size,D))
	mol_com=np.sum(mxyz,axis=1)/(mol_mass)

	# translate com of molecules for plus delta and minus delta shift
	#create a nmolecules x mol_size x D array
	xyz=xyz.reshape((-1,mol_size,D))

	xyz_p=mol_trans(xyz,delta,mol_com)
	xyz_m=mol_trans(xyz,-delta,mol_com)
	
	
	#write out shifted molecules to lammps dump files
	write_lammps_dump(xyz_p,lammps_dump,lammps_dump_header,xyz_p_file)
	write_lammps_dump(xyz_m,lammps_dump,lammps_dump_header,xyz_m_file)
	
	
	#find interface of surface slab
	
	#read in coordinates of surface from lammps dump file
	xyz_surface,mass_surface,lammps_dump_surface=read_lammps_dump('xyz_surface.txt')
	
	#wrap surface coordinates to box
	xyz_surface=wrapcoord(xyz_surface,boxD,box_center=[0,0,0])
	
	#calculate distribution of surfaces as a function of z
	surface_freq,bins=np.histogram(xyz_surface[:,2],range(-int(boxD[2]/2.0),int(boxD[2]/2.0)+1,1))
	
	#find min and max of surface surface, cutoff used arbitrarily to avoid pure vapor
	#assumes that there is a value of more than 20 surfaces in an 80x80x1 A slab
	vapor_cut=.003125*boxD[0]*boxD[1]
	minl=np.argmax(surface_freq>vapor_cut)
	minr=np.size(surface_freq)-np.argmax(surface_freq[::-1]>vapor_cut)
	
	#estimate bulk surface number density from middle 20 A of surface slab
	slabcenter=(minl+minr)/2.0
	surface_density=np.average(surface_freq[slabcenter-10:slabcenter+10])
	
	#find interface
	interface_l=np.argmax(surface_freq>surface_density/2.0)
	interface_r=np.size(surface_freq)-np.argmax(surface_freq[::-1]>surface_density/2.0)-1
	
	#find the interface of the monomer and air
	#wrap monomer molecules 
	mol_com_wrapped=wrapcoord(mol_com,boxD,box_center=[0,0,0])
	mol_shift=mol_com_wrapped-mol_com

	xyz_wrap=np.zeros_like(xyz)
	for i in range(np.size(xyz,axis=0)):
		xyz_wrap[i]=xyz[i]+mol_shift[i]
		
	#find monomers who have first headgroup within cutoff of surface surface
	#find z coord of interfaces, this is estimated as the middle of the bin
	z_l=np.mean(bins[interface_l:interface_l+2])
	z_r=np.mean(bins[interface_r:interface_r+2])
	
	#returns indicies of monomer molecules that are within cutoff of surface surface
	monomer_l=np.where(abs(z_l-xyz_wrap[:,0,2])<surface_cut)[0]
	monomer_r=np.where(abs(z_r-xyz_wrap[:,0,2])<surface_cut)[0]

	monomer_list=np.concatenate((monomer_l,r)) #HERE
	#remove any duplicates. There should not be any duplicates anyway.
	list=np.unique(monomer_list)
	xyz_monomer=xyz_wrap[list]
	
	#compute distribution of non hydophilic monomer beads
	freq,bins=np.histogram(xyz_monomer[:,2:,2].flatten(),range(-int(boxD[2]/2.0),int(boxD[2]/2.0)+1,1))

	#Estimate bulk monomer bead density
	#0.00125 is probably too low for production
	vapor_cut=0.00125*boxD[0]*boxD[1]
	minl1_monomer=np.argmax(freq>monomer_vapor_cut)
	minl2_monomer= - (np.argmax(freq[:int(slabcenter)][::-1]>monomer_vapor_cut)-np.size(monomer_freq[:int(slabcenter)]))-1
	minlcenter=int((minl1_monomer+minl2_monomer)/2.0)
	
	minr2_monomer=np.size(freq)-np.argmax(monomer_freq[::-1]>monomer_vapor_cut)-1
	minr1_monomer=np.argmax(freq[int(slabcenter):]>monomer_vapor_cut)+np.size(monomer_freq)-np.size(monomer_freq[:int(slabcenter)])-1
	minrcenter=int((minr1_monomer+minr2_monomer)/2.0)
	
	density=np.average(monomer_freq[minlcenter-3:minlcenter+3])
	density+=np.average(monomer_freq[minrcenter-3:minrcenter+3])
	density=monomer_density/2.0
	
	
	interface_l=np.argmax(monomer_freq>monomer_density/2.0)
	interface_r=np.size(monomer_freq)-np.argmax(monomer_freq[::-1]>monomer_density/2.0)-1
	#number of monomer beads between vapor and surface interface
	count=np.sum(monomer_freq[interface_monomer_l:interface_l+1])+np.sum(monomer_freq[interface_r:interface_monomer_r+1])
	#volume of monomer 
	vol=monomer_count*(4.0/3.0)*np.pi*3.7244**3.0 

	#volume of interface
	#assuemes bins all have same width as first bin
	interface_vol=(abs(bins[interface_l]-bins[interface_l])+abs(bins[interface_monomer_r]-bins[interface_r]))*abs(bins[1]-bins[0])*boxD[0]*boxD[1]
	
	f_volfrac=open('temp_vol_frac.txt','w')
	f_volfrac.write(str(vol/interface_vol)+' '+str(interface_vol)+'\n\n')
	f_volfrac.close()
	return 0

a=main()

