#!/usr/bin/env python  

import os  
import Symmetric_Monomer_Addition 
import sys

f1 = open(sys.argv[4],'a') #Open in read format                                                          
out_stream = f1

lt_file = sys.argv[2]

os.system('src/packmol/packmol < monomer_additional.inp')

os.system('src/moltemp_src.2013-3-03/moltemplate.sh  -nocheck -xyz monomer_additional.xyz %s'%(lt_file))

os.system('rm -rf output_ttree')
os.system('rm system_trimmed.i*')
os.system('rm system_conf.i*')

#Writing an output file 
p = "-------------------------------------------------------------------"
p += "\nlower interface: %5.3f; upper interface: %5.3f\n"  %(Symmetric_Monomer_Addition.DefineSurface.lower_interface ,Symmetric_Monomer_Addition.DefineSurface.upper_interface)
p+= "top layer has %i monomers and bottom layer has %i monomers\n" %(Symmetric_Monomer_Addition.DefineSurface.total_top, Symmetric_Monomer_Addition.DefineSurface.total_bottom)
p+= "-------------------------------------------------------------------\n\n"

out_stream.write(p)
out_stream.close()

