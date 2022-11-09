# To generate topology for cellulose/matrices
# Combined both cellulose and matrix in one code.
# Users can choose Martini V2.0 or Martini V3.0 for cellulose
# Version 2.0: Sept-19-2022
# Version 3.0: Oct-31-2022
# I/p cmd (polymer matrix): python gen_top.py matrix_name nmons
# matrix_name - pla, p3hb, petg 
# I/p cmd (cellulose): python gen_top.py gro_fname cell_dp ncnf_per_bundle nbundle
# acetfrac martini_ver molname
# Authors: Shalini Jayaraman Rukmani & Vaidyanathan Sethuraman
import os
import sys
import numpy
import re
import shutil
import glob
import math
import subprocess
import auxgen_polytop
import auxgen_martinibeads
import auxgen_celltop_v2
import auxgen_celltop_v3
#----Read input file - matrix/filename,ncnf_per_bundle,acetfrac-----------
beads_per_mon = 4
if len(sys.argv) == 3: #polymer matrix - design and exit
    print('Generating CG polymer for ', sys.argv[1])
    if sys.argv[1] == 'petg':
        auxgen_polytop.design_petg(int(sys.argv[2]),sys.argv[1].upper())
    elif sys.argv[1] == 'pla':
        auxgen_polytop.design_pla(int(sys.argv[2]),sys.argv[1].upper())
    elif sys.argv[1] == 'p3hb':
        auxgen_polytop.design_p3hb(int(sys.argv[2]),sys.argv[1].upper())
    else:
        raise RuntimeError('Unknown matrix input: ' + sys.argv[1])
    exit()
elif len(sys.argv) != 8: # Cellulose - read parameters
    print('Requires exactly 8 arguments instead of ', len(sys.argv),\
          str(sys.argv))
    raise RuntimeError('Unknown number of arguments')
#---------------- Analyzing cellulose -----------------------------
print('Input file name for cellulose CG: ',sys.argv[1])

#-----Process input data - Cellulose/Acetylated Cellulose-----------
fname        = str(sys.argv[1])   # Coarse-grained file
cell_dp      = int(sys.argv[2])   # degree of polymerization of cellulose
ch_per_cnf   = int(sys.argv[3])   # Number of cnf chains per bundle
ncnf_bundles = int(sys.argv[4])   # num of cellulose bundles
acetfrac     = float(sys.argv[5]) # acetylated fraction
martini_ver  = int(sys.argv[6])   # Martini version
molname      = str(sys.argv[7])   # Molecular name
#------Check Martini version-------------------------------------------
if martini_ver not in [2,3]:
    raise RuntimeError("Unknown Martini version: " + str(martini_ver))

#-----Check input file is present-----------------------------------
if not os.path.exists(fname):
    raise RuntimeError(fname + ' not found in path!')

#-----Generate log file---------------------------------------------
fout = auxgen_martinibeads.gen_logfile(fname,ncnf_bundles,acetfrac,\
                                       cell_dp,ch_per_cnf)
residarr,resnamearr,aidarr,anamearr,rxarr,\
    ryarr,rzarr,massarr = auxgen_martinibeads.read_gro_file(fname)

if acetfrac == 0: # Sanity check is not possible for acetylated systems
    if len(rxarr) != ncnf_bundles*ch_per_cnf*cell_dp*beads_per_mon:
        raise RuntimeError('Unequal number of atoms betweem inputs and gro file: '\
                           + str(ncnf_bundles*ch_per_cnf*cell_dp*beads_per_mon)\
                           + '\t' + str(len(rxarr)))
#-----Generate bead list--------------------------------------------
glycan_list = auxgen_martinibeads.create_martini_beads(cell_dp,ncnf_bundles,\
                                                       ch_per_cnf,residarr,\
                                                       aidarr,anamearr)
#-----Generate topology----------------------------------------------
if martini_ver == 2:
    # Create bond list
    bond_list = auxgen_celltop_v2.create_bond_list(cell_dp,ncnf_bundles,ch_per_cnf,glycan_list)
    # Create angle list
    angle_list = auxgen_celltop_v2.create_angle_list(cell_dp,ncnf_bundles,ch_per_cnf,glycan_list)
    # Create dihedral list
    dihed_list = auxgen_celltop_v2.create_dihedral_list(cell_dp,ncnf_bundles,ch_per_cnf,glycan_list)
    # COMBINE AND WRITE
    posre_fname = str(fname.split("/")[-1]).split(".gro")[0] 
    ftop = auxgen_celltop_v2.write_celltop(posre_fname,residarr,\
                                           resnamearr,aidarr,anamearr,\
                                           massarr,bond_list,angle_list,dihed_list) 

elif martini_ver == 3: # Parameters for cellulose are added in cell_martini3.itp
    # COMBINE AND WRITE
    posre_fname = str(fname.split("/")[-1]).split(".gro")[0] 
    auxgen_celltop_v3.GLCB14(cell_dp,ncnf_bundles,ch_per_cnf,glycan_list,\
                             posre_fname,molname)
#-----End code----------------------------------------------------
