# To generate psf files 
# Main file: gen_top.py
# Reference: Martini3-Polyply
import os
import sys
import numpy
import re
import shutil
import glob
import math
import subprocess
#----------------------------------------------------------------------------------

def renumber_CEL_beads(fbnd, file_ctr, atom_list_for_lengths, bond_list, new_bond_list):
   for j in range(0, len(bond_list)):
      if file_ctr >  1:
         val1 = bond_list[j][0] + len(atom_list_for_lengths)
         val2 = bond_list[j][1] + len(atom_list_for_lengths)
      else:
         val1 = bond_list[j][0]
         val2 = bond_list[j][1]
      fbnd.write('%d\t%d\n' %(val1,val2))
      new_bond_list.append([val1,val2])
   file_ctr += 1
   return new_bond_list
#----------------------------------------------------------------------------------
def renumber_poly_beads(fbnd, nbeads_per_chain, nchains, polymer_start_index, bo_list, new_bond_list):
    for chain_count in range(0, nchains):
        for j in range(0, len(bo_list)):
            val1 = bo_list[j][0] + polymer_start_index + chain_count*(nbeads_per_chain) 
            val2 = bo_list[j][1] + polymer_start_index + chain_count*(nbeads_per_chain)
            fbnd.write('%d\t%d\n' %(val1, val2))
            new_bond_list.append([val1,val2])
        fbnd.write('next chain')
    return new_bond_list
#----------------------------------------------------------------------------------
### Write PSF file
def write_psf(fname,residarr,resnamearr,aidarr,anamearr,massarr,new_bond_list):
    segnamearr = []
    ftop = open(fname,'w')
    charge = 0.0 # Martini polymers and cellulose have zero charge
    rem_str = "; qtot 0"
    ftop.write("PSF CMAP\n\n")
    ftop.write("     3 !NTITLE\n") # three lines after this line
    ftop.write("REMARKS generated structure x-plor psf file from auxgen_top.py\n")
    ftop.write("REMARKS based on topology files: CELLULOSE_acetylated.top and <polymer>.itp\n")
    ftop.write("REMARKS <polymer> can pla, pp, p3hb and petg\n\n")
    natoms = len(anamearr)
    ftop.write("    %d !NATOM\n" %(natoms))
    last_column = 0
    seg_ctr = 1
    segname = 'U' + str(seg_ctr)
    segnamearr.append(segname)
    for i in range(1, len(aidarr)):
        if residarr[i] >= residarr[i-1]:
            segname = 'U' + str(seg_ctr)
            segnamearr.append(segname)
        else:
            seg_ctr += 1
            segname = 'U' + str(seg_ctr)
            segnamearr.append(segname)

    for i in range(0, len(aidarr)):
        ftop.write("%8s %4s %4s %4s %4s %4s %14s%14s%8s\n" \
                   %(str(aidarr[i]),segnamearr[i],str(residarr[i]),\
                     resnamearr[i],anamearr[i],anamearr[i],str(charge),\
                     str(massarr[i]),str(last_column)))
    ftop.write("\n")
    nbonds = len(new_bond_list)
    ftop.write("   %d !NBOND: bonds\n" %(nbonds))
    for i in range(0, len(new_bond_list)):
        ftop.write("%10s %10s" %(str(new_bond_list[i][0]),str(new_bond_list[i][1])))
        if i%4 == 3:
            ftop.write("\n")
    ftop.write("\n")
#----------------------------------------------------------------------------------
# if __name__
if __name__ == '__main__':
    main()
#----------------------------------------------------------------------------------
