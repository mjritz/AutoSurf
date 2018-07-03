#! /usr/bin/env/ python

import sys
import numpy as np
import pdb_create
import variables

#Writes bond and angle information for nonpolar (NP), polar (P) and terminal polar (TP) groups

structure = open(sys.argv[1], 'r')
in_file  = open('monomer.xyz', 'r')
out_file = open('monomer_struct.lt', 'wa')

RDF, xbox, ybox, zbox, surface_filename, surface_mol, pack_tol, bead_dist, boxx, boxy, boxz,     bead_per, main_branch, equil_time_short,equil_time_long, monomer_x, monomer_y, monomer_zlow, monomer_zhi,         per_loop, steps,temp, z_distance, epsilon, sigma, r_distance, equil_dcd,        prod_dcd, cutoff = variables.read_variables('input/input_file.txt')


#Initialize string that will be appended to file with bond information
s = []
ilines, i,j,k,l = 0,0,0,0,0
ID = []
xyz = []
y = []
z = []

#Write beginning of file 
s = ("""                                                                                                 
import "mie_ff.lt"                                                                                   

Monomer {                                                                                                   
                                                                                                             
  write('Data Atoms') {                                                                                  
""")   

for ilines in in_file.readlines():  
  if ilines.startswith ('NP') or ilines.startswith('P') or ilines.startswith('TP'): 
    data = ilines.split()
    if data[0] == 'NP':
      s+="    $atom:NP%i $mol:. @atom:MieMonomer/%s  %3.2f  %8.4f  %8.4f  %8.4f\n" %(i+1,data[0],float(0.0),float(data[1]), float(data[2]),float(data[3]))
      i+=1
    if data[0] == 'P':
      s+="    $atom:P%i $mol:. @atom:MieMonomer/%s  %2.2f  %8.4f  %8.4f %8.4f\n"%(j+1,data[0],float(0.0),float(data[1]), float(data[2]),float(data[3]))
      j+=1
    if data[0] == 'TP':
      s+="    $atom:TP%i $mol:. @atom:MieMonomer/%s  %2.2f  %8.4f  %8.4f %8.4f\n"%(k+1,data[0],float(0.0),float(data[1]), float(data[2]),float(data[3]))
      k+=1
s+="  }"       

in_file  = open('monomer.xyz', 'r') 

#Initialize string that will be appended to file with bond information                                   
head_groups = 0                                                                                          
number_of_atoms = 0                                                                                      
sum_side_beads = 0                                                                                       
branches = []                                                                                               
branch_location =[]                                                                                      
g = []                                                                                                   
e,n,i,j,k,l,m = 0,0,0,0,0,0,0                                                                                
                                                                                                         
g ="""                                                                                                   
                                                                                                         
  write ('Data Bonds') {                                                                                 
"""          

for ilines in structure.readlines():                                                                     
  if 'B' in ilines or 'T' in ilines or 'H' in ilines or 'E' in ilines:                                                       
    atoms_on_each_line = len(ilines.split())                                                             
    number_of_atoms += atoms_on_each_line                                                                
  if ilines.count('B') >= 0:                                                                             
    number_side_beads = len(ilines.split())-1                                                            
    branches.append(number_side_beads)                                                                   
    sum_side_beads += number_side_beads                                                                  
  if ilines.count('H') or ilines.count('E') > 0:                                                                              
    head_group_beads = len(ilines.split())                                                               
    head_groups += head_group_beads                                                                      
main_tail = number_of_atoms-sum_side_beads-head_groups   


for x in range(main_tail):                                                                               
  if k<main_tail-1:                                                                                      
    g += "    $bond:bond%i @bond:MieMonomer/%s%s $atom:%s%i $atom:%s%i\n"%(k+1,'NP','NP','NP',n+1,'NP',n+2)
  if k == main_tail-1:                                                                                   
    g+= "    $bond:bond%i @bond:MieMonomer/%s%s $atom:%s%i $atom:%s%i\n"%(k+1,'NP','P','NP',n+1,'P',i+1)
  n+=1                                                                                                   
  k += 1 

for x in range(head_groups-1):               
  if k<main_tail-1+head_groups-1:                                                                                    
    g += "    $bond:bond%i @bond:MieMonomer/%s%s $atom:%s%i $atom:%s%i\n"%(k+1,'P','P','P',i+1,'P',i+2)
  if k == main_tail-1+head_groups-1:
    g+= "    $bond:bond%i @bond:MieMonomer/%s%s $atom:%s%i $atom:%s%i\n"%(k+1,'P','TP','P',i+1,'TP',e+1)
    e+1
  i+=1
  k+=1

for lines in in_file.readlines():                                                                      
  if lines.startswith ('NP') or lines.startswith('P'):                                                
    data = lines.split()                                                                           
    if data[2] !=  '0.00000':                                                                      
      branch_location.append(data[2])    

for x in range(len(branches)):
  for x in range(branches[l]):                                                                           
    if branch_location[m] == '%3.5f' %(float(bead_dist)):                                                                  
      g+= "    $bond:bond%i @bond:MieMonomer/%s%s $atom:%s%i $atom:%s%i\n"%(k+1,'NP','NP','NP',l+1,'NP',n+1)
      n-=1                                                                                               
    elif branch_location[m] == '-%3.5f' %(float(bead_dist)):                                                               
      g+= "    $bond:bond%i @bond:MieMonomer/%s%s $atom:%s%i $atom:%s%i\n"%(k+1,'NP','NP','NP',n+1,'NP',l+1) 
    else:                                                                                                
      g+= "    $bond:bond%i @bond:MieMonomer/%s%s $atom:%s%i $atom:%s%i\n"%(k+1,'NP','NP','NP',n+1,'NP',n+2)
    n+=1                                                                                                 
    m+=1                                                                                                 
    k+=1                                                                                                 
  if branches[l] != 0:                                                                                   
    n+=1                                                                                                 
  l+=1                                                                                                   
g+="  }" 
g+="""

  write ('Data Angles'){
"""

n,i,j,k,l,m,o = 0,0,0,0,0,0,0                                                                                

for x in range(main_tail):                                                                               
  if k<main_tail-2:                                                                                      
    g += "    $angle:angle%i @angle:MieMonomer/%s%s%s $atom:%s%i $atom:%s%i $atom:%s%i\n"%(k+1,'NP','NP','NP','NP',n+1,'NP',n+2,'NP',n+3)
  n+=1                                                                                                   
  k += 1 

for lines in in_file.readlines():                                                                      
  if lines.startswith ('NP') or lines.startswith('P'):                                                
    data = lines.split()                                                                           
    if data[2] !=  '0.00000':                                                                      
      branch_location.append(data[2])    
k-=2

for x in range(len(branches)):
  for x in range(branches[l]):                                                                           
    if branch_location[m] == '%3.5f' %(float(bead_dist)):
      if l < main_tail-2: 
        g+= "    $angle:angle%i @angle:MieMonomer/%s%s%s $atom:%s%i $atom:%s%i $atom:%s%i\n"%(k+1,'NP','NP','NP','NP',n+1,'NP',l+1,'NP',l+2)
        k+=1
      if l > 1:
        g+= "    $angle:angle%i @angle:MieMonomer/%s%s%s $atom:%s%i $atom:%s%i $atom:%s%i\n"%(k+1,'NP','NP','NP','NP',n+1,'NP',l+1,'NP',l)
        k+=1
    elif branch_location[m] == '-%3.5f' %(float(bead_dist)):                                                               
      if l < main_tail-2:
        g+= "    $angle:angle%i @angle:MieMonomer/%s%s%s $atom:%s%i $atom:%s%i $atom:%s%i\n"%(k+1,'NP','NP','NP','NP',n+1,'NP',l+1,'NP',l+2)
        k+=1
      if l> 1:
        g+= "    $angle:angle%i @angle:MieMonomer/%s%s%s $atom:%s%i $atom:%s%i $atom:%s%i\n"%(k+1,'NP','NP','NP','NP',n+1,'NP',l+1,'NP',l)
        k+=1
    elif branch_location[m] == '-%3.5f' %(2*float(bead_dist)):
        g+= "    $angle:angle%i @angle:MieMonomer/%s%s%s $atom:%s%i $atom:%s%i $atom:%s%i\n"%(k+1,'NP','NP','NP','NP',n+1,'NP',n+2,'NP',l+1)
        k+=1
    elif branch_location[m] == '%3.5f' %(2*float(bead_dist)):
        g+= "    $angle:angle%i @angle:MieMonomer/%s%s%s $atom:%s%i $atom:%s%i $atom:%s%i\n"%(k+1,'NP','NP','NP','NP',n+1,'NP',n,'NP',l+1)
        k+=1
    else:                                                                                                
      if branch_location[m] < '-%3.5f' %(2*float(bead_dist)):
        g+= "    $angle:angle%i @angle:MieMonomer/%s%s%s $atom:%s%i $atom:%s%i $atom:%s%i\n"%(k+1,'NP','NP','NP','NP',n+1,'NP',n+2,'NP',n+3)
        k+=1
      if branch_location[m] > '%3.5f' %(2*float(bead_dist)):
        g+= "    $angle:angle%i @angle:MieMonomer/%s%s%s $atom:%s%i $atom:%s%i $atom:%s%i\n"%(k+1,'NP','NP','NP','NP',n+1,'NP',n+2,'NP',n+3)
        k+=1
    n+=1                                                                                                 
    m+=1                                                                                                 
  l+=1                                                                                                   

g+= """ }                                                                                                
                                                                                                                     
}"""       
    
out_file.write(s)
out_file.write(g)
out_file.close()

if sum_side_beads == 0:
    descript = 'linear'
if sum_side_beads != 0:
    descript = 'branched'
