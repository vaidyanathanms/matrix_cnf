 # To generate psf file for cellulose/polymer matrices composite
# Version 1.0: Oct-2-2022
# I/p cmd: python gen_psf.py matrix_name nmons nchains fname_list ncnf acetfrac cell_dp ch_per_cnf fname
# matrix_name - pla, p3hb, petg 
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
from auxgen_psf import *
#from auxgen_celltop_v2 import *
#from auxgen_celltop_v3 import *

#----Read input file - matrix/filename,ncnf_fibers,acetfrac-----------
ch_per_cnf   = 18 # default
if len(sys.argv) >= 5 : 
    nmons = int(sys.argv[2])
    nchains = int(sys.argv[3])
    matrix_name = str(sys.argv[1])
    print('Generating CG polymer for ', sys.argv[1])
    if sys.argv[1] == 'petg':
        at_list,bo_list,an_list,di_list = design_petg(int(sys.argv[2]),sys.argv[1].upper())
    elif sys.argv[1] == 'pla':
        at_list,bo_list,an_list,di_list = design_pla(int(sys.argv[2]),sys.argv[1].upper())
    elif sys.argv[1] == 'p3hb':
        at_list,bo_list,an_list,di_list = design_p3hb(int(sys.argv[2]),sys.argv[1].upper())
    else:
        raise RuntimeError('Unknown matrix input: ' + sys.argv[1])
else:
    print("incorrect number of arguments")
    exit()
if len(sys.argv) == 10:
    ch_per_cnf = int(sys.argv[8]) # number of chains per bundle

# Open file to write bond list
fbnd = open('bondlist.txt','w')
#---------------- Analyzing cellulose -----------------------------
if len(sys.argv) > 5 :
   print('Input file name for cellulose CG: ',sys.argv[4])
   #-----Process input data - Cellulose/Acetylated Cellulose-----------
   fname    = str(sys.argv[4]) # Coarse-grained file
   ncnf     = int(sys.argv[5]) # num of cellulose chain in a bundle
   acetfrac = float(sys.argv[6]) # acetylated fraction
   cell_dp  = int(sys.argv[7]) # degree of polymerization of cellulose

   #-----Check input file is present-----------------------------------
   fname = sorted(glob.glob(str(fname)))
   bond_list = []
   new_bond_list = []
   atom_list = []
   atom_list_for_lengths = []
   segidarr=[]
   file_ctr = 1 # to later check for the first file
   natoms = 0
   nmol = 0
   for f in fname:
       if not os.path.exists(f):
           raise RuntimeError(f + ' not found in path!')
   
       #-----Generate log file---------------------------------------------
       #fout = gen_logfile(f,ncnf,acetfrac,cell_dp,ch_per_cnf)
       residarr,resnamearr,aidarr,anamearr,rxarr,ryarr,rzarr,massarr = read_gro_file(f) # read only cellulose gro file
       #-----Generate bead list--------------------------------------------
       glycan_list = create_martini_beads(cell_dp,ncnf,ch_per_cnf,residarr,\
                                          aidarr,anamearr)
       # Append atom list
       atom_list.append(aidarr)
       # get atom list length for each CNF
       if file_ctr >  1:
          for elem in atom_list[file_ctr-2]:
              atom_list_for_lengths.append(int(elem))
       # Create bond list
       bond_list = create_bond_list(cell_dp,ncnf,ch_per_cnf,glycan_list)
       # For renumbering beads
       new_bond_list = renumber_CEL_beads(fbnd, file_ctr, atom_list_for_lengths, bond_list, new_bond_list)     
       file_ctr += 1
   for i in range(0, len(atom_list)):
       natoms += len(atom_list[i]) 
   polymer_start_index = natoms 
   new_bond_list = renumber_poly_beads(fbnd, len(at_list), nchains, polymer_start_index, bo_list, new_bond_list)
else:
   new_bond_list = []
   polymer_start_index = 0
   new_bond_list = renumber_poly_beads(fbnd, len(at_list), nchains, polymer_start_index, bo_list, new_bond_list)

#------Generate PSF file--------------------------------------
if len(sys.argv) > 5:
   fname=str(sys.argv[9])
   outfile_name = "CNF-polymer.psf"
else:
   fname=str(sys.argv[4])
   outfile_name = "polymer.psf" 
residarr,resnamearr,aidarr,anamearr,rxarr,ryarr,rzarr,massarr = read_gro_file(fname) # read cellulose-polymer gro file
fpsf = write_psf(outfile_name,residarr,resnamearr,aidarr,anamearr,massarr,new_bond_list)
