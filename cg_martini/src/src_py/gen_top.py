# To generate topology for cellulose/matrices
# Combined both cellulose and matrix in one code.
# Users can choose Martini V2.0 or Martini V3.0 for cellulose
# Version 2.0: Sept-19-2022
# I/p cmd (polymer matrix): python gen_top.py matrix_name nmons
# matrix_name - pla, p3hb, petg 
# I/p cmd (cellulose): python gen_top.py fname ncnf acetfrac cell_dp martini_ver ch_per_cnf
# ch_per_cnf is optional (default is 18)
# Authors: Shalini Jayaraman Rukmani & Vaidyanathan Sethuraman
import os
import sys
import numpy
import re
import shutil
import glob
import math
import subprocess
from auxgen_polytop import *
from auxgen_celltop_v2 import *
from auxgen_celltop_v3 import *

#----Read input file - matrix/filename,ncnf_fibers,acetfrac-----------
ch_per_cnf   = 18 # default
martini_ver  = -1 # default
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
elif len(sys.argv) == 7: # Cellulose - read parameters
    martini_ver = int(sys.argv[5]) # Martini version
    ch_per_cnf  = int(sys.argv[6]) # number of chains per bundle
elif len(sys.argv) != 6: 
    print('Unknown number of arguments: ', len(sys.argv),\
          str(sys.argv))
    exit()

if martini_ver not in [-1,2,3]:
    raise RuntimeError("Unknown Martini version: " + str(martini_ver))
#---------------- Analyzing cellulose -----------------------------
print('Input file name for cellulose CG: ',sys.argv[1])

#-----Process input data - Cellulose/Acetylated Cellulose-----------
fname    = str(sys.argv[1]) # Coarse-grained file
ncnf     = int(sys.argv[2]) # num of cellulose bundles
acetfrac = float(sys.argv[3]) # acetylated fraction
cell_dp  = int(sys.argv[4]) # degree of polymerization of cellulose

#-----Check input file is present-----------------------------------
if not os.path.exists(fname):
    raise RuntimeError(fname + ' not found in path!')

#-----Generate log file---------------------------------------------
fout = gen_logfile(fname,ncnf,acetfrac,cell_dp,ch_per_cnf)
residarr,resnamearr,aidarr,anamearr,rxarr,ryarr,rzarr,massarr = read_gro_file(fname)
                                                                
#-----Generate bead list--------------------------------------------
if martini_ver == 2:
    glycan_list = auxgen_celltop_v2.create_martini_beads(cell_dp,ncnf,ch_per_cnf,residarr,\
                                       aidarr,anamearr)
    # Create bond list
    bond_list = auxgen_celltop_v2.create_bond_list(cell_dp,ncnf,ch_per_cnf,glycan_list)
    # Create angle list
    angle_list = auxgen_celltop_v2.create_angle_list(cell_dp,ncnf,ch_per_cnf,glycan_list)
    # Create dihedral list
    dihed_list = auxgen_celltop_v2.create_dihedral_list(cell_dp,ncnf,ch_per_cnf,glycan_list)
    # COMBINE AND WRITE
    posre_fname = str(fname.split("/")[-1]).split(".gro")[0] 
    ftop = auxgen_celltop_v2.write_celltop(posre_fname,residarr,\
                                           resnamearr,aidarr,anamearr,\
                                           massarr,bond_list,angle_list,dihed_list) 

elif martini_ver == 3: # Parameters for cellulose are added in cell_martini3.itp
    glycan_list = auxgen_celltop_v3.create_martini_beads(cell_dp,ncnf,ch_per_cnf,residarr,\
                                       aidarr,anamearr)
    # Create bond list
    bond_list = auxgen_celltop_v3.create_bond_list(cell_dp,ncnf,ch_per_cnf,glycan_list)
    # Create angle list
    angle_list = auxgen_celltop_v3.create_angle_list(cell_dp,ncnf,ch_per_cnf,glycan_list)
    # Create dihedral list
    dihed_list = auxgen_celltop_v3.create_dihedral_list(cell_dp,ncnf,ch_per_cnf,glycan_list)
    # COMBINE AND WRITE
    posre_fname = str(fname.split("/")[-1]).split(".gro")[0] 
    ftop = auxgen_celltop_v3.write_celltop(posre_fname,residarr,\
                                           resnamearr,aidarr,anamearr,\
                                           massarr,bond_list,angle_list,dihed_list) 
#-----End code----------------------------------------------------
