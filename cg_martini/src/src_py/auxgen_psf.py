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
def read_gro_file(inpfyle):
    fmt_at = b'5s5s5s5s8s8s8s'
    aidarr = []; anamearr = [];residarr = []; resnamearr = []
    rxarr = []; ryarr = []; rzarr = []; massarr = []
    with open(inpfyle,'r') as fin:
        fin.readline()
        nbeads = int(fin.readline().strip())
        for line in fin:
            if line.startswith(';'): # skip lines starting with #
                continue
            if not line: # skip empty lines
                continue
            elif len(line.strip().split()) == 3 or len(line.strip().split()) == 9:
                lbox = line.strip().split()
            else:
                residarr.append(int(line[:5].strip()))
                resnamearr.append(line[5:10].strip())
                anamearr.append(line[10:15].strip())
                aidarr.append(int(line[15:20].strip()))
                rxarr.append(float(line[20:28].strip()))
                ryarr.append(float(line[28:36].strip()))
                rzarr.append(float(line[36:44].strip()))

            if 'S' == line[10]: # small
                mval = 54
            elif 'T' == line[10]: # tiny
                mval = 36
            else: # Default Martini regular
                mval = 72
            massarr.append(mval)
    if len(residarr) == 0:
        raise RuntimeError('No atoms read - file corrupted?')
    return residarr,resnamearr,aidarr,anamearr,rxarr,ryarr,rzarr,massarr
#----------------------------------------------------------------------------------
# Creating MARTINI beads list
def create_martini_beads(cell_dp,ncnf,ch_per_cnf,residarr,\
                         aidarr,anamearr):
    # Create empty array
    glycan_list = [[] for i in range(cell_dp*ncnf*ch_per_cnf)]
    # Sort beads
    headid_ptr = -1; resid_gro = 0; glycan_cnt = -1
    while resid_gro < len(residarr):
        if residarr[resid_gro] == headid_ptr: #check if resID is same
            while residarr[resid_gro] == headid_ptr:
                glycan_list[glycan_cnt].append(aidarr[resid_gro])
                glycan_list[glycan_cnt].append(anamearr[resid_gro])
                resid_gro += 1 #move to next "atom" in gro file
                if resid_gro == len(residarr): #break while loop
                    break
        else: #update new glycan
            glycan_cnt += 1 #move to new "glycan"
            headid_ptr = residarr[resid_gro] #update pointer to new glycan
    return glycan_list

#----------------------------------------------------------------------------------
# Create bond list between the following glycan type names
# Beta(1-4)
#      1          7
#     /          / \  
#   4  -  2 -  8  -  6
#     \ /        \
#      3          5 
def create_bond_list(cell_dp,ncnf,ch_per_cnf,glycan_list):
    # Create empty list
    bond_list = [[]]
    bnd_cntr = 0
    fbnd = open('bond.txt','w')
    for i in range(len(glycan_list)):
        for j in range(0,len(glycan_list[i]),2):
            for k in range(j+2,len(glycan_list[i]),2):
                if (glycan_list[i][j+1] == 'T1' and glycan_list[i][k+1] == 'S4'):
                   bond_list[bnd_cntr].append(glycan_list[i][j])
                   bond_list[bnd_cntr].append(glycan_list[i][k])
                   bnd_cntr += 1
                   fbnd.write('%d\t%d\n' %(glycan_list[i][j],glycan_list[i][k]))
                   bond_list.append([])
                elif (glycan_list[i][j+1] == 'T2' and glycan_list[i][k+1] == 'R3'):
                   bond_list[bnd_cntr].append(glycan_list[i][j])
                   bond_list[bnd_cntr].append(glycan_list[i][k])
                   bnd_cntr += 1
                   fbnd.write('%d\t%d\n' %(glycan_list[i][j],glycan_list[i][k]))
                   bond_list.append([])
                elif (glycan_list[i][j+1] == 'R3' and glycan_list[i][k+1] == 'S4'): 
                   bond_list[bnd_cntr].append(glycan_list[i][j])
                   bond_list[bnd_cntr].append(glycan_list[i][k])
                   bnd_cntr += 1
                   fbnd.write('%d\t%d\n' %(glycan_list[i][j],glycan_list[i][k]))
                   bond_list.append([])
                elif (glycan_list[i][j+1] == 'T2' and glycan_list[i][k+1] == 'S4'):
                   bond_list[bnd_cntr].append(glycan_list[i][j])
                   bond_list[bnd_cntr].append(glycan_list[i][k])
                   bnd_cntr += 1
                   fbnd.write('%d\t%d\n' %(glycan_list[i][j],glycan_list[i][k]))
                   bond_list.append([])

    # Connection between glycans (T2-S4, T2-T2 and S4-S4)
    for ch in range(ncnf*ch_per_cnf):
        for glyres in range(cell_dp-1):
            i = cell_dp*ch + glyres
            j = glycan_list[i].index('T2') - 1; # -1 to point to ID of PX
            k = glycan_list[i+1].index('S4') - 1
            bond_list[bnd_cntr].append(glycan_list[i][j])
            bond_list[bnd_cntr].append(glycan_list[i+1][k])
            bnd_cntr += 1
            fbnd.write('%d\t%d\n' %(glycan_list[i][j],glycan_list[i+1][k]))
            bond_list.append([])
    for ch in range(ncnf*ch_per_cnf):
        for glyres in range(cell_dp-1):
            i = cell_dp*ch + glyres
            j = glycan_list[i].index('S4') - 1; # -1 to point to ID of the corresponding bead
            k = glycan_list[i+1].index('S4') - 1
            bond_list[bnd_cntr].append(glycan_list[i][j])
            bond_list[bnd_cntr].append(glycan_list[i+1][k])
            bnd_cntr += 1
            fbnd.write('%d\t%d\n' %(glycan_list[i][j],glycan_list[i+1][k]))
            bond_list.append([])
    for ch in range(ncnf*ch_per_cnf):
        for glyres in range(cell_dp-1):
            i = cell_dp*ch + glyres
            j = glycan_list[i].index('T2') - 1; # -1 to point to ID of the corresponding bead
            k = glycan_list[i+1].index('T2') - 1
            bond_list[bnd_cntr].append(glycan_list[i][j])
            bond_list[bnd_cntr].append(glycan_list[i+1][k])
            bnd_cntr += 1
            fbnd.write('%d\t%d\n' %(glycan_list[i][j],glycan_list[i+1][k]))
            bond_list.append([])
    del(bond_list[-1])
    return bond_list

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
    print (nchains, nbeads_per_chain, len(bo_list), polymer_start_index)
    fbnd.write('start+poly')
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
