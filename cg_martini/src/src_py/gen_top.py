# To generate topology for cellulose/matrices

import os
import sys
import numpy
import re
import shutil
import glob
import math
import subprocess
from auxgen_top import *

#----Read input file - filename, ncnf_fibers, acetfrac--------------
ch_per_cnf   = 18 # default
if len(sys.argv) == 6:
    ch_per_cnf = int(sys.argv[5]) # number of chains per bundle
elif len(sys.argv) != 5: 
    print('Unknown number of arguments: ', len(sys.argv),\
          str(sys.argv))
    exit()

print('Input file name: ',sys.argv[1])

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


# CREATE BOND LIST

# CREATE ANGLE LIST
# CREATE DIHEDRAL LIST
# CREATE ATOMTYPE LIST
# CREATE NONBONDED PARAMETER LIST
# CREATE BONDED PARAMETER LIST
# COMBINE AND WRITE

