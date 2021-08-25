# To check guassian polymer matrix chains
# Add above to PACKMOL along with modified cellulose chains or additional
# polymers

import os
import sys
import numpy
import re
import shutil
import glob
import math
import subprocess
from aux_pack import * # function definitions
#------------------------------------------------------------------

# Input data
matrix   = 'pla' #pla/pp
mat_fyle = 'pla_inp.pdb' # matrix input file - ONLY PDB
cnf_fyle = 'modifiedCNF_m1' # cellulose/mod cellulose input
nmons    = 65 # number of matrix monomers
nchains  = 50 # number of matrix chains
ncnf     = 1 # number of cellulose bundles
mod_cell = 1 # 1 - modified cellulose; 0 - native cellulose
add_poly = 'None'
inppack  = 'pack_cellulose.inp' # PACKMOL input file
gaus_tol = 0.05 # tolerance for checking gaussianity
fin_box  = 1.5 # final box size relative to max dimension of cnf/matrix
run_pack = 0 # 1-run packmol
#------------------------------------------------------------------

# Directory paths
main_dir  = os.getcwd() # current dir
cnf_dir   = '/home/v0e/allcodes/files_cnf/make_acetylated_cellulose'
pack_dir  = '/home/v0e/packmol'
if not os.path.isdir(cnf_dir): # cellulose directory
    print("FATAL ERROR: ", cnf_dir, " not found")
    exit("Check cnf directory path")
scr_dir   = '/lustre/or-scratch/cades-bsd/v0e' # scratch dir
if not os.path.isdir(scr_dir):
    print("FATAL ERROR: ", scr_dir, " not found\n")
    exit("Check scratch directory path\n")
if not os.path.isdir(pack_dir) and run_pack == 1:
    print("FATAL ERROR: ", pack_dir, "not found\n")
    exit("Check packmol directory path\n")
    if not os.path.exists(pack_dir + '/packmol'):
        print("FATAL ERROR: Packmol executable not found\n")
        exit("Check packmol executable\n")
scr_dir  = scr_dir + '/cellulose_nanofibers'
if fin_box < 1:
    raise RuntimeError("Unphysical final box magnification\n")
#------------------------------------------------------------------

# Main code

# Check for polymer matrix files and compute their dimensions
check_inp_files(main_dir,mat_fyle) #polymer matrix
matdir,rgmax = check_gaussianity_and_write(mat_fyle,nmons,nchains,matrix,gaus_tol)

# Check for cnf dimensions and compare with rg of polymer matrix
xmin,ymin,zmin,xmax,ymax,zmax = find_cnf_dim(cnf_dir,main_dir,cnf_fyle)
dmax = max(xmax-xmin,ymax-ymin,zmax-zmin,rgmax)

# Start writing PACKMOL files
fpack   = open(inppack,'w')
packout = packmol_headers(fpack,matrix)

# Pack CNF/matrix/additional polymers
pack_cellulose_chains(fpack,cnf_fyle,ncnf,xmin,ymin,zmin,xmax,ymax,zmax,fin_box,dmax)
pack_polymer_matrix(matdir,matrix,nchains,xmin,ymin,zmin,xmax,ymax,zmax,fpack,fin_box,dmax)
#if add_poly != 'None':
#    pack_extra_polymers()

# Close and run PACKMOL
fpack.close()
if run_pack:
    run_packmol(inppack,pack_dir)

# Create final directories in scratch
make_fin_dirs_and_copy(scr_dir,matrix,mod_cell,run_pack,packout)

    



