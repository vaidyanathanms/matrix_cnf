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

# Input directory paths
main_dir = os.getcwd() # current dir
nativ_cnf= '/home/v0e/allcodes/files_cnf/elementary_fibrils' #native cellulose dir
acet_dir = '/home/v0e/allcodes/files_cnf/make_acetylated_cellulose' #acet_cnf dir
pack_exe = '/home/v0e/packmol/packmol' # packmol executable
poly_mat = '/home/v0e/allcodes/files_cnf/polymer_matrices' #i/o dir poly matrices
scr_dir  = '/lustre/or-scratch/cades-bsd/v0e' # scratch dir

# Input data - Polymer matrix
matrix   = 'pla' #pla/pp/petg/p3hb
mat_pdb  = 'step3_input.pdb' # matrix input pdb file - ONLY PDB
nmons    = 40 # number of matrix monomers
nchains  = 81 # number of matrix chains
gaus_tol = 0.05 # tolerance for checking gaussianity

# Input data - Cellulose/Acetylated Cellulose/Additives
acet_val = 1 # m1 - 1, m3 - 3, m7 - 7, m11 - 11
acet_per = 0.5 # fraction of acetylated cellulose
acet_new = 1 # 0 - use old, 1-delete and regenerate
ncnf     = 1 # number of cellulose bundles (18 chains)
cell_dp  = 20 # degree of polymerization of cellulose chains
acet_tol = 0.1 # tolerance for acetylation
acetpref = 'modified_m' # prefix for acetylated files
acet_att = 100 # maximum attempts to create acetylated cellulose
add_poly = 'None'

# Input data - Packmol
inppack  = 'pack_cellulose.inp' # PACKMOL input file
fin_box  = 1.1 # final box size relative to max dimension of cnf/matrix
run_pack = 1 # 1-run packmol
packsh   = 'run_packmol_pyinp.sh'
#------------------------------------------------------------------

# Check for directory paths and input consistency
pack_sup = poly_mat + '/cnf_packed_' + matrix + '_N_' + str(nmons)\
           + '_M_' + str(nchains)  # final packed output super dir
poly_dir = poly_mat + '/' + matrix + '/charmm_' + matrix + '_N_' + str(nmons)\
             + '_M_' + str(nchains) # input polymer matrix dir
gmx_mat  = poly_dir + '/gromacs'

check_dir(nativ_cnf)
check_dir(acet_dir)
check_dir(poly_mat)
check_dir(poly_dir)
check_dir(gmx_mat)
check_dir(scr_dir)

if fin_box < 1:
    raise RuntimeError("Unphysical final box magnification\n")
#------------------------------------------------------------------

# Make output directories/files
scr_dir   = scr_dir + '/cellulose_nanofibers'
if acet_per != 0:
    acet_fyle = acetpref + str(acet_val)
    mod_cell  = 1
else:
    acet_fyle = native_cell_dp + str(cell_dp)
    mod_cell = 0
#------------------------------------------------------------------

# Main code

print('Begin creating packmol files...')

print('Begin checking input files and create output directories...')

# Check for polymer matrix files 
check_inp_files(gmx_mat,mat_pdb) 

# Create output directories
pack_mat = create_output_dirs(pack_sup,acet_val,acet_per,add_poly)

print('Checking matrix input files for gaussian chains...')
# Check for Gaussian input chains
matdir,rgmax = check_gaussianity_and_write(gmx_mat,mat_pdb,nmons,\
                                           nchains,matrix,gaus_tol,pack_mat)

# Acetylate chains if needed
if acet_per != 0:
    print('Making acetylated chains ..')
    make_acet_cell(acet_dir,acet_val,cell_dp,acet_per,acet_tol,\
                   acet_fyle,acet_att,pack_mat,acet_new,nativ_cnf)
    os.chdir(main_dir)

print('Determining optimal dimensions of the box...')
# Check for cnf dimensions and compare with rg of polymer matrix
xmin,ymin,zmin,xmax,ymax,zmax = find_cnf_dim(acet_dir,acet_fyle,pack_mat)
dmax = max(xmax-xmin,ymax-ymin,zmax-zmin,rgmax)

# Start writing PACKMOL files
packfyle = pack_mat + '/' + inppack
fpack   = open(packfyle,'w')
packout = packmol_headers(fpack,matrix,pack_mat,acet_per,acet_val,add_poly)

# Pack CNF/matrix/additional polymers
print('Writing packmol scripts for packing cellulose and matrix chains...')
pack_cellulose_chains(fpack,pack_mat,acet_fyle,ncnf,xmin,ymin,zmin,xmax,ymax,zmax,fin_box,dmax)
pack_polymer_matrix(matdir,matrix,nchains,xmin,ymin,zmin,xmax,ymax,zmax,fpack,fin_box,dmax)
#if add_poly != 'None':
#    pack_extra_polymers()

fpack.close() # close PACKMOL input file

# Run PACKMOL
if run_pack:
    run_packmol(packfyle,pack_exe,pack_mat,packsh,main_dir)

# Combine psf/top files for the system
#combine_psf_top_files()

# Clean up psf/pdb files
#clean_and_sort_files()

# Generate shell input files
#generate_sh_files()

os.chdir(main_dir)
# Create final directories in scratch
# print('Creating input directories of MD simulations...')
# make_fin_dirs_and_copy(scr_dir,matrix,mod_cell,run_pack,packout)
    
print('Done :) ..')



